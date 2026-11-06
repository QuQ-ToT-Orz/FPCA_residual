#### 0 Packages ####
# Robustness of the FPC2-mortality association, with the multilevel fPCA adjustments added: does day-to-day variation (level 2)
# carry information beyond each participant's usual shape (level 1)?
# Level-2 scores are per day;
# L2_SD1 to L2_SD3: the SD over the seven days of the three leading level-2 scores, and
# L2_total overall day-to-day variability: sqrt(sum over all level-2 components of the score variance)
gc()
library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(survey)

#### 1 Settings ####
dir_path <- "./results/mortality/"
load(paste0(dir_path, "cohorts_weighted.RData"))
base_terms <- c("Age", "Gender", "Race", "Income", "MESOR_z", "Amplitude_z", "cosinor_phase_sin", "cosinor_phase_cos")
fpc_terms <- c("FPC1", "FPC2", "FPC3")
main_terms <- c(base_terms, fpc_terms)
activity_terms <- c("wear_z", "asinh_TAC", "MVPA_z")
outcomes <- c(mortstat = "All-cause", cvd_death = "Cardiovascular", cancer_death = "Cancer")
methods <- c("unweighted", "weighted")
cohort_label <- c(count = "Hip counts, 2003-2006", mims = "Wrist MIMS, 2011-2014", pooled = "Pooled (fixed effect)")

# MFPCA level-2 summaries per participant
z <- function(x) as.numeric(scale(x))
cosinor_fit <- "weekly"
l2_sd_terms <- c("L2_SD1", "L2_SD2", "L2_SD3")
level2_describe <- list()
for (ds in names(cohorts)) {
  suffix <- if (ds == "mims") "new" else "old"
  load(paste0("./data/", ds, "/Act_Analysis_asinh_", suffix, ".RData"))
  day_id <- as.character(arrange(Act_Analysis_asinh, SEQN, WEEKDAY)$SEQN)
  rm(Act_Analysis_asinh)
  load(paste0("./results/cosinor_residual_fpca/", ds, "/residual_fpca_", cosinor_fit, "_", suffix, ".RData"))
  S2 <- as.matrix(mfpca_fit$scores$level2)
  # SD over the seven days of each level-2 score; 
  # summed score variances are the integrated day-to-day variance of the level-2 curve.
  sds <- apply(S2, 2, function(s) tapply(s, day_id, sd))
  level2 <- data.frame(SEQN = rownames(sds), L2_SD1 = sds[, 1], L2_SD2 = sds[, 2], L2_SD3 = sds[, 3], L2_total = sqrt(rowSums(sds^2)))
  ev1 <- mfpca_fit$evalues$level1
  ev2 <- mfpca_fit$evalues$level2
  level2_describe[[ds]] <- data.frame(dataset = ds, level2_components = ncol(S2), level1_share = sum(ev1) / (sum(ev1) + sum(ev2)), level2_comp1_3_share = sum(ev2[1:3]) / sum(ev2))
  # Per cohort SD
  cohorts[[ds]]$d <- cohorts[[ds]]$d %>%
    left_join(level2, by = "SEQN") %>%
    mutate(across(c(all_of(l2_sd_terms), "L2_total"), z))
  rm(subject_residual, fitted_subject_residual, fpca_fit, fpca_variance, mfpca_fit, mfpca_variance, S2)
  gc()
}
level2_describe <- bind_rows(level2_describe)
print(level2_describe, row.names = FALSE)
# Correlation of FPC2 with the level-2 summaries
for (ds in names(cohorts)) {
  cat(ds, ": Spearman rho with FPC2\n")
  print(round(cor(cohorts[[ds]]$d$FPC2, cohorts[[ds]]$d[, c(l2_sd_terms, "L2_total", "acro_csd")], method = "spearman"), 2))
}

# Correlation of the Model 1 variables with the level-2 summaries
corr_vars <- c(setdiff(main_terms, "Race"), "IS", "IV", "acro_csd", l2_sd_terms, "L2_total")
corr <- lapply(cohorts, function(x) {
  d <- x$d[corr_vars]
  d$Income[d$Income == "Missing"] <- NA
  d[] <- lapply(d, as.numeric)
  cor(d, method = "spearman", use = "pairwise.complete.obs")
})
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
         dataset = factor(cohort_label[dataset], cohort_label[names(corr)]))
p_corr <- ggplot(corr_long, aes(var2, var1, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(abs(rho) >= 0.7, sprintf("%.2f", rho), "")), size = 2, colour = "grey15") +
  facet_wrap(~dataset) +
  scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", midpoint = 0, limits = c(-1, 1),
                       name = "Spearman rho") +
  coord_fixed() +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(x = NULL, y = NULL)
p_corr

ext_terms <- c("ext_minimum", "ext_alpha", "ext_log_beta")   # Extended cosinor beyond the ordinary cosinor
rhythm_terms <- c("IS", "IV", "acro_csd")   # acro_csd: circular SD of the daily acrophases
lifestyle_terms <- c("BMI", "SmokeCigs", "DrinkStatus", "EducationAdult", "MobilityProblem")
combined_terms <- c(ext_terms, rhythm_terms, "MobilityProblem")
specs <- list(
  list(group = "Adjustment", label = "Model 1", terms = main_terms),
  list(group = "Adjustment", label = "+ extended cosinor", terms = c(main_terms, ext_terms)),
  list(group = "Adjustment", label = "+ extended cosinor, IS, IV and acrophase SD", terms = c(main_terms, ext_terms, rhythm_terms)),
  list(group = "Adjustment", label = "+ mobility", terms = c(main_terms, "MobilityProblem")),
  list(group = "Adjustment", label = "+ all above (Model 2)", terms = c(main_terms, combined_terms)),
  list(group = "Adjustment", label = "Model 2 + wear time", terms = c(main_terms, combined_terms, "wear_z")),
  # list(group = "Adjustment", label = "+ total activity", terms = c(main_terms, "asinh_TAC")),
  list(group = "Adjustment", label = "Model 2 + MVPA", terms = c(main_terms, combined_terms, "MVPA_z")),
  list(group = "Adjustment", label = "Model 2 + BMI, smoking, alcohol, education",terms = c(main_terms, combined_terms, lifestyle_terms)),
  list(group = "Adjustment", label = "Model 2 + all above", terms = c(main_terms, combined_terms, lifestyle_terms, "wear_z", "MVPA_z")),
  list(group = "Day-to-day", label = "Model 2 + level-2 SD1-SD3", terms = c(main_terms, combined_terms, l2_sd_terms)),
  list(group = "Day-to-day", label = "Model 2 + overall day-to-day variability", terms = c(main_terms, combined_terms, "L2_total")),
  list(group = "Sensitivity", label = "Deaths in first 2 years excluded", terms = main_terms, early = TRUE))
all_terms <- unique(unlist(lapply(specs, `[[`, "terms")))
level2_terms <- c(l2_sd_terms, "L2_total")

#### 2 Model functions ####
designs <- lapply(cohorts, function(x) svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~w_ipw,
                                                 nest = TRUE, data = x$d))
cox <- function(d, terms, event, method, ds) {
  f <- as.formula(paste0("Surv(permth_exm / 12, ", event, ") ~ ", paste(terms, collapse = " + ")))
  if (method == "unweighted") return(coxph(f, data = d, ties = "efron"))
  svycoxph(f, design = subset(designs[[ds]], SEQN %in% d$SEQN))
}
estimate <- function(f, term) {
  b <- unname(coef(f)[term])
  se <- sqrt(unname(vcov(f)[term, term]))
  c(b = b, se = se, p = 2 * pnorm(-abs(b / se)))
}

#### 3 Fit every specification ####
# FPC2 in every model; the level-2 terms of the day-to-day models are kept too
robust <- list()
for (ev in names(outcomes)) for (ds in names(cohorts)) {
  d <- cohorts[[ds]]$d
  d <- d[complete.cases(d[, c("permth_exm", ev, all_terms)]), ]
  for (m in methods) {
    for (s in specs) {
      dd <- d
      if (isTRUE(s$early)) dd <- dd[!(dd$mortstat == 1 & dd$permth_exm < 24), ]
      f <- cox(dd, s$terms, ev, m, ds)
      for (tm in c("FPC2", intersect(s$terms, level2_terms))) {
        k <- estimate(f, tm)
        robust[[length(robust) + 1]] <- data.frame(
          outcome = ev, dataset = ds, method = m, group = s$group, label = s$label, term = tm,
          n = nrow(dd), deaths = sum(dd[[ev]]),
          b = k[["b"]], se = k[["se"]], p = k[["p"]],
          ph_p = if (m == "unweighted") cox.zph(f)$table[tm, "p"] else NA_real_) # proportional-hazards diagnostics
      }
    }
  }
}
robust <- bind_rows(robust)

# Pooled across cohorts: fixed-effect inverse-variance; heterogeneity: Cochran's Q
pooled <- robust %>%
  group_by(outcome, method, group, label, term) %>%
  summarise(het_p = pchisq(sum((b - sum(b / se^2) / sum(1 / se^2))^2 / se^2), 1, lower.tail = FALSE), # Heterogeneity
            b = sum(b / se^2) / sum(1 / se^2), se = sqrt(1 / sum(1 / se^2)),
            n = sum(n), deaths = sum(deaths), .groups = "drop") %>%
  mutate(dataset = "pooled", p = 2 * pnorm(-abs(b / se)))
robust <- bind_rows(robust, pooled) %>% mutate(hr = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se))

# Level-2 terms themselves (HR per cohort SD of within-person variability)
level2_hr <- robust %>%
  filter(term %in% level2_terms) %>%
  mutate(hr_ci = sprintf("%.2f (%.2f-%.2f)", hr, lo, hi), p = round(p, 3)) %>%
  select(method, outcome, label, term, dataset, hr_ci, p) %>%
  pivot_wider(names_from = dataset, values_from = c(hr_ci, p))
print(as.data.frame(level2_hr), row.names = FALSE)

#### 4 Forest plots ####
robust_plot <- function(m) {
  r <- robust %>%
    filter(method == m, term == "FPC2") %>%
    mutate(label = factor(label, rev(unique(sapply(specs, `[[`, "label")))),
           group = factor(group, unique(sapply(specs, `[[`, "group"))),
           outcome = factor(outcomes[outcome], outcomes),
           cohort = factor(cohort_label[dataset], cohort_label),
           sig = ifelse(p < 0.05, "p < 0.05", "p >= 0.05"))
  # position_dodge stacks the first level at the bottom; reverse so hip is on top and pooled at the bottom
  ggplot(r, aes(hr, label, colour = cohort, shape = cohort, alpha = sig, group = forcats::fct_rev(cohort))) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.6, position = position_dodge(width = 0.75)) +
    geom_point(size = 2.2, position = position_dodge(width = 0.75)) +
    facet_grid(group ~ outcome, scales = "free_y", space = "free_y") +
    scale_x_log10(breaks = c(0.5, 0.8, 1, 1.25, 1.5, 2, 3)) +
    scale_colour_manual(values = c("#2360A0", "#B74F27", "grey15")) +
    scale_shape_manual(values = c(16, 16, 18)) +
    scale_alpha_manual(values = c(1, 0.35), guide = "none") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", strip.text.y = element_text(angle = 270)) +
    labs(title = paste(if (m == "weighted") "Survey-weighted" else "Unweighted", "Cox models"),
         x = "Hazard ratio per SD of FPC2 (log scale)", y = NULL, colour = NULL, shape = NULL)
}
p_weighted <- robust_plot("weighted")
p_unweighted <- robust_plot("unweighted")
p_weighted
p_unweighted
