# Participant characteristics, 
# the leading residual components, 
# FPC2 interpretation (groups and correlations) 
#### 0 Packages ####
gc()
library(dplyr)
library(tidyr)
library(survey)

#### 1 Settings ####
cosinor_fit <- "weekly"
dir_path <- "./results/mortality/"
load(paste0(dir_path, "cohorts_weighted.RData"))
base_terms <- c("Age", "Gender", "Race", "Income", "MESOR_z", "Amplitude_z", "cosinor_phase_sin", "cosinor_phase_cos")
fpc_terms <- c("FPC1", "FPC2", "FPC3")
nf <- function(x) x / sqrt(sum(x^2))
clock <- function(minutes) sprintf("%02d:%02d", floor(minutes / 60) %% 24, round(minutes %% 60) %% 60)
mean_sd <- function(x, digits = 1) sprintf(paste0("%.", digits, "f (%.", digits, "f)"), mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
pct <- function(x) sprintf("%.1f", 100 * mean(x, na.rm = TRUE))
med_iqr <- function(x, digits = 1) {
  q <- quantile(x, c(0.5, 0.25, 0.75), na.rm = TRUE)
  sprintf(paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"), q[1], q[2], q[3])
}
clock_iqr <- function(x) {
  q <- quantile(x, c(0.5, 0.25, 0.75), na.rm = TRUE)
  sprintf("%s (%s-%s)", clock(q[1]), clock(q[2]), clock(q[3]))
}

#### 2 Participant characteristics ####
characteristics <- function(d) c(
  "Participants" = as.character(nrow(d)),
  "Age, years" = mean_sd(d$Age),
  "Female, %" = pct(d$Gender == "Female"),
  "Race: White, %" = pct(d$Race == "White"),
  "Race: Black, %" = pct(d$Race == "Black"),
  "Race: Hispanic, %" = pct(d$Race == "Hispanic"),
  "Race: other, %" = pct(d$Race == "Other"),
  "Income-to-poverty ratio < 1.3" = pct(d$Income[d$Income != "Missing"] == "Low"),
  "Income-to-poverty ratio > 3.5, %" = pct(d$Income[d$Income != "Missing"] == "High"),
  "BMI, kg/m2" = mean_sd(d$BMI),
  "Current smoker, %" = pct(d$SmokeCigs == "Current"),
  "More than high school, %" = pct(d$EducationAdult == "More than high school"),
  "Any mobility difficulty, %" = pct(d$MobilityProblem == "Any Difficulty"),
  "Mean activity per minute" = mean_sd(d$TAC / 1440),
  "MVPA, median (IQR) min/day" = med_iqr(d$MVPA, 0),
  "Wear time, hours/day" = mean_sd(d$WT / 60),
  "Cosinor MESOR" = mean_sd(d$cosinor_mes, 2),
  "Cosinor amplitude" = mean_sd(d$cosinor_amp, 2),
  "Cosinor acrophase, median (IQR)" = clock_iqr(d$cosinor_acrotime),
  "Interdaily stability (IS)" = mean_sd(d$IS, 2),
  "Intradaily variability (IV)" = mean_sd(d$IV, 2),
  "Follow-up, median (IQR) years" = med_iqr(d$permth_exm / 12),
  "Deaths, n" = as.character(sum(d$mortstat)),
  "Cardiovascular deaths, n" = as.character(sum(d$cvd_death)),
  "Cancer deaths, n" = as.character(sum(d$cancer_death)))
table_demo <- data.frame(characteristic = names(characteristics(cohorts$count$d)),
                     hip = characteristics(cohorts$count$d), wrist = characteristics(cohorts$mims$d), row.names = NULL)
table_demo

#### 3 Leading residual components ####
# Shape: share of each eigenfunction explained by a sinusoid with period 24, 12 or 8 hours. 
wave_share <- function(phi, period) {
  X <- apply(cbind(cos(2 * pi * (1:1440) / period), sin(2 * pi * (1:1440) / period)), 2, nf)
  sum(crossprod(X, nf(phi))^2)
}
# Similarity: inner product of the hip and wrist eigenfunctions
similarity <- sapply(1:3, function(k) sum(nf(cohorts$count$efunctions[, k]) * nf(cohorts$mims$efunctions[, k])))
components <- list()
for (ds in names(cohorts)) {
  suffix <- if (ds == "mims") "new" else "old"
  e <- new.env()
  load(paste0("./results/cosinor_residual_fpca/", ds, "/residual_fpca_", cosinor_fit, "_", suffix, ".RData"), envir = e)
  share <- e$fpca_variance
  phi <- cohorts[[ds]]$efunctions
  components[[ds]] <- data.frame(dataset = ds, component = fpc_terms, variance_share = share[1:3],
  cumulative_share = cumsum(share)[1:3],
  period_24h = apply(phi, 2, wave_share, 1440),
  period_12h = apply(phi, 2, wave_share, 720),
  period_8h = apply(phi, 2, wave_share, 480),
  hip_wrist_similarity = similarity,
  n_components_90 = which(cumsum(share) >= 0.9)[1])
}
components <- bind_rows(components)
components %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
components

#### 4 FPC2 interpretation ####
# Participants with FPC2 <= -1 SD (low), between (middle) and >= +1 SD (high)
fpc2_groups <- bind_rows(lapply(names(cohorts), function(ds) {
  cohorts[[ds]]$d %>%
    mutate(group = factor(ifelse(FPC2 >= 1, "High", ifelse(FPC2 <= -1, "Low", "Middle")), c("Low", "Middle", "High"))) %>%
    group_by(group) %>%
    summarise(n = n(), age = mean(Age), female_pct = 100 * mean(Gender == "Female"),
    mobility_pct = 100 * mean(MobilityProblem == "Any Difficulty"),
    activity_per_min = mean(TAC / 1440), mvpa_median = median(MVPA),
    mesor = mean(cosinor_mes), amplitude = mean(cosinor_amp),acrophase_median = clock(median(cosinor_acrotime)), IS = mean(IS), IV = mean(IV), died_pct = 100 * mean(mortstat), .groups = "drop") %>%
    mutate(dataset = ds, .before = 1)
}))
fpc2_groups %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))

#### 5 Save ####
out_dir <- paste0(dir_path, "interpretation/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(table_demo, paste0(out_dir, "participant_characteristics.csv"), row.names = FALSE)
write.csv(components, paste0(out_dir, "residual_components.csv"), row.names = FALSE)
write.csv(fpc2_groups, paste0(out_dir, "fpc2_groups.csv"), row.names = FALSE)