#### 0 Packages ####
gc()
library(dplyr)
library(haven)
library(ggplot2)

#### 1 Settings ####
cosinor_fit <- "weekly"
z <- function(x) as.numeric(scale(x))
nf <- function(x) x / sqrt(sum(x^2))
expanded <- TRUE   # FALSE for the default run
out_suffix <- if (expanded) "_expanded" else ""
dir_path <- paste0("./results/mortality", out_suffix, "/")
dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)

#### 2 Build one cohort ####
build_cohort <- function(dataset) {
  suffix <- paste0(if (dataset == "mims") "new" else "old", out_suffix)
  cycles <- if (dataset == "mims") c("G", "H") else c("C", "D")

  load(paste0("./data/", dataset, "/data_analysis_", suffix, ".RData"))
  data_analysis$SEQN <- as.character(data_analysis$SEQN)
  load(paste0("./results/cosinor_residual_fpca/", dataset, "/residual_fpca_", cosinor_fit, "_", suffix, ".RData"))
  subject_id <- rownames(fpca_fit$scores)

  # Weekly cosinor and extended cosinor parameters
  load(paste0("./results/cosinor_analysis/", dataset, "/weekly_cosinor_parameters_", suffix, ".RData"))
  cosinor <- weekly_parameters$cosinor %>%
    transmute(SEQN = as.character(SEQN), cosinor_mes = mes, cosinor_amp = amp,
    cosinor_acrotime = acrotime,
    cosinor_phase_sin = sin(2 * pi * acrotime / 1440), cosinor_phase_cos = cos(2 * pi * acrotime / 1440))
  extended <- weekly_parameters$extended %>%
    transmute(SEQN = as.character(SEQN), ext_minimum = minimum, ext_amp = amp, ext_alpha = alpha,
    ext_log_beta = log(beta),
    ext_phase_sin = sin(2 * pi * acrotime / 1440), ext_phase_cos = cos(2 * pi * acrotime / 1440))

  load(paste0("./data/", dataset, "/Act_Analysis_asinh_", suffix, ".RData"))
  volume <- Act_Analysis_asinh %>%
    mutate(SEQN = as.character(SEQN)) %>%
    group_by(SEQN) %>%
    summarise(TAC = mean(TAC), MVPA = mean(MVPA), WT = mean(WT), .groups = "drop")
  rm(Act_Analysis_asinh)

  # Family income-to-poverty ratio
  income <- bind_rows(lapply(paste0("DEMO_", cycles), function(f) {
    read_xpt(paste0("./data/nhanes_demo/", f, ".xpt"))[, c("SEQN", "INDFMPIR")]})) %>%
    transmute(SEQN = as.character(as.integer(SEQN)), INDFMPIR, Income = factor(case_when(
                is.na(INDFMPIR) ~ "Missing", 
                INDFMPIR < 1.3 ~ "Low",
                INDFMPIR < 3.5 ~ "Middle", 
                TRUE ~ "High"), c("Low", "Middle", "High", "Missing")))

  # Top three components
  scores <- data.frame(SEQN = subject_id, fpca_fit$scores[, 1:3])
  names(scores)[-1] <- paste0("FPC", 1:3)
  
  d <- data.frame(SEQN = subject_id) %>%
    left_join(data_analysis, by = "SEQN") %>%
    left_join(cosinor, by = "SEQN") %>%
    left_join(extended, by = "SEQN") %>%
    left_join(volume, by = "SEQN") %>%
    left_join(income, by = "SEQN") %>%
    left_join(scores, by = "SEQN")

  # Cause-specific deaths (NCHS leading cause), for cause-specific Cox models:
  # cardiovascular = 1 (heart disease) or 5 (cerebrovascular disease); 
  # cancer = 2 (malignant neoplasms).
  cause <- as.integer(as.character(d$ucod_leading))
  d <- d %>% mutate(
    cvd_death = as.integer(mortstat == 1 & cause %in% c(1, 5)),
    cancer_death = as.integer(mortstat == 1 & cause %in% 2),
    # Per cohort SD
    MESOR_z = z(cosinor_mes), Amplitude_z = z(cosinor_amp),
    asinh_TAC = z(asinh(TAC)), MVPA_z = z(MVPA), wear_z = z(WT),
    across(paste0("FPC", 1:3), z))

  # Eigenfunctions are kept to align component signs between cohorts
  list(d = d, efunctions = fpca_fit$efunctions[, 1:3])
}

#### 3 Build both cohorts and align component signs ####
cohorts <- list(count = build_cohort("count"), mims = build_cohort("mims"))
align <- function(efun, prefix) {
  s <- sapply(1:3, function(k) sum(nf(cohorts$count[[efun]][, k]) * nf(cohorts$mims[[efun]][, k])))
  for (k in which(s < 0)) {
    cohorts$mims$d[[paste0(prefix, k)]] <<- -cohorts$mims$d[[paste0(prefix, k)]]
    cohorts$mims[[efun]][, k] <<- -cohorts$mims[[efun]][, k]
  }
  abs(s)
}
similarity <- align("efunctions", "FPC")
cat("Hip vs wrist eigenfunction similarity, FPC1-FPC3:", round(similarity, 3), "\n")
# Hip vs wrist eigenfunction similarity, FPC1-FPC3: 0.95 0.959 0.681
# (expanded) Hip vs wrist eigenfunction similarity, FPC1-FPC3: 0.949 0.946 0.696 
for (ds in names(cohorts)) {
  d <- cohorts[[ds]]$d
  cat(ds, ": n =", nrow(d), "| deaths =", sum(d$mortstat), "| CVD deaths =", sum(d$cvd_death), "| cancer deaths =", sum(d$cancer_death), "\n")
}
# count : n = 2615 | deaths = 463 | CVD deaths = 112 | cancer deaths = 101 
# mims : n = 3208 | deaths = 312 | CVD deaths = 106 | cancer deaths = 77

# expanded
# count : n = 2917 | deaths = 507 | CVD deaths = 120 | cancer deaths = 109 
# mims : n = 3605 | deaths = 362 | CVD deaths = 125 | cancer deaths = 91 

#### 4 Spearman Correlation ####
corr_vars <- c("Age", "MESOR_z", "Amplitude_z", 
               "cosinor_phase_sin", "cosinor_phase_cos", 
               "ext_minimum", "ext_amp", "ext_alpha", "ext_log_beta", "ext_phase_sin", "ext_phase_cos",
               paste0("FPC", 1:3),
               "IS", "IV", "asinh_TAC", "MVPA_z", "wear_z")
corr <- lapply(cohorts, function(x) cor(x$d[corr_vars], method = "spearman", use = "pairwise.complete.obs"))
high_corr <- bind_rows(lapply(names(corr), function(ds) {
  r <- corr[[ds]]
  i <- which(upper.tri(r) & abs(r) >= 0.7, arr.ind = TRUE)
  data.frame(dataset = ds, var1 = rownames(r)[i[, 1]], var2 = colnames(r)[i[, 2]], rho = r[i])
})) %>% arrange(dataset, desc(abs(rho)))
print(high_corr %>% mutate(rho = round(rho, 2)), row.names = FALSE)

corr_long <- bind_rows(lapply(names(corr), function(ds) {
  r <- corr[[ds]]
  i <- which(lower.tri(r), arr.ind = TRUE)
  data.frame(dataset = ds, var1 = rownames(r)[i[, 1]], var2 = colnames(r)[i[, 2]], rho = r[i])
})) %>%
  mutate(var1 = factor(var1, rev(corr_vars)), var2 = factor(var2, corr_vars),
         dataset = factor(c(count = "Hip counts, 2003-2006", mims = "Wrist MIMS, 2011-2014")[dataset]))
ggplot(corr_long, aes(var2, var1, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(abs(rho) >= 0.7, sprintf("%.2f", rho), "")), size = 1.5, colour = "grey15") +
  facet_wrap(~dataset) +
  scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", midpoint = 0, limits = c(-1, 1),
                       name = "Spearman rho") +
  coord_fixed() +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(x = NULL, y = NULL)

save(cohorts, file = paste0(dir_path, "cohorts.RData"))
