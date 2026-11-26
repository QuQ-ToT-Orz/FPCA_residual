#### Comparison with standard FPCA ####
# FPC2 and its standard counterpart (FPC4) in the mortality models #
# (a) Masking: each score entered alone with age, sex, ethnicity and income, without and with the cosinor parameters.
# (b) Robustness: the adjustment sets. standard FPC1-FPC3 cannot be entered with the cosinor parameters
gc()
library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(survey)

#### 1 Settings ####
dir_path <- "./results/mortality/"
load(paste0(dir_path, "cohorts_weighted.RData"))
base_terms <- c("Age", "Gender", "Race", "Income", "MESOR_z", "Amplitude_z", "cosinor_phase_sin", "cosinor_phase_cos", "wear_z")
fpc_terms <- c("FPC1", "FPC2", "FPC3")
main_terms <- c(base_terms, fpc_terms)
activity_terms <- c("asinh_TAC", "MVPA_z")
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
  nf <- function(x) x / sqrt(sum(x^2))
  phi1 <- mfpca_fit$efunctions$level1
  l1_similarity <- sapply(1:3, function(k) abs(sum(nf(phi1[, k]) * nf(fpca_fit$efunctions[, k]))))
  l1_score_cor <- sapply(1:3, function(k) {
    s <- sign(sum(phi1[, k] * cohorts[[ds]]$efunctions[, k]))
    cor(s * mfpca_fit$scores$level1[cohorts[[ds]]$d$SEQN, k], cohorts[[ds]]$d[[paste0("FPC", k)]])
  })
  level2_describe[[ds]] <- data.frame(dataset = ds, level2_components = ncol(S2), level1_share = sum(ev1) / (sum(ev1) + sum(ev2)), level2_comp1_3_share = sum(ev2[1:3]) / sum(ev2),
    similarity_FPC1 = l1_similarity[1], similarity_FPC2 = l1_similarity[2], similarity_FPC3 = l1_similarity[3],
    score_cor_FPC1 = l1_score_cor[1], score_cor_FPC2 = l1_score_cor[2], score_cor_FPC3 = l1_score_cor[3])
  # Per cohort SD
  cohorts[[ds]]$d <- cohorts[[ds]]$d %>%
    left_join(level2, by = "SEQN") %>%
    mutate(across(c(all_of(l2_sd_terms), "L2_total"), z))
  rm(subject_residual, fitted_subject_residual, fpca_fit, fpca_variance, mfpca_fit, mfpca_variance, S2)
  gc()
}

ext_terms <- c("ext_minimum", "ext_alpha", "ext_log_beta")   # Extended cosinor beyond the ordinary cosinor
rhythm_terms <- c("IS", "IV", "acro_csd")   # acro_csd: circular SD of the daily acrophases
lifestyle_terms <- c("BMI", "SmokeCigs", "DrinkStatus", "EducationAdult", "MobilityProblem")
combined_terms <- c(ext_terms, rhythm_terms, "MobilityProblem")
disease_terms <- c("Diabetes", "CHF", "CHD", "Stroke",  "Cancer")   # Self-reported at baseline (DIQ010, MCQ160B/C/F, MCQ220)
clinical_terms <- c("SYS", "TC_HDL")   # Systolic BP, total/HDL cholesterol ratio
health_terms <- c(disease_terms, clinical_terms)
specs <- list(
  list(group = "Adjustment", label = "Model 1", terms = main_terms),
  list(group = "Adjustment", label = "+ extended cosinor", terms = c(main_terms, ext_terms)),
  list(group = "Adjustment", label = "+ extended cosinor, IS, IV and acrophase SD", terms = c(main_terms, ext_terms, rhythm_terms)),
  list(group = "Adjustment", label = "+ mobility", terms = c(main_terms, "MobilityProblem")),
  list(group = "Adjustment", label = "+ all above (Model 2)", terms = c(main_terms, combined_terms)),
  # list(group = "Adjustment", label = "+ total activity", terms = c(main_terms, "asinh_TAC")),
  list(group = "Adjustment", label = "Model 2 + MVPA", terms = c(main_terms, combined_terms, "MVPA_z")),
  list(group = "Adjustment", label = "Model 2 + BMI, smoking, alcohol, education",terms = c(main_terms, combined_terms, lifestyle_terms)),
  list(group = "Adjustment", label = "Model 2 + all above", terms = c(main_terms, combined_terms, lifestyle_terms, "MVPA_z")),
  list(group = "Health status", label = "Model 1 +\nCHF, CHD, stroke, diabetes, cancer +,\nSBP, total/HDL cholesterol ratio", terms = c(main_terms, health_terms)),
  list(group = "Health status", label = "Model 2 +\nCHF, CHD, stroke, diabetes, cancer +,\nSBP, total/HDL cholesterol ratio", terms = c(main_terms, combined_terms, health_terms)),
  list(group = "Day-to-day", label = "Model 2 + level-2 SD1-SD3", terms = c(main_terms, combined_terms, l2_sd_terms)),
  list(group = "Day-to-day", label = "Model 2 + overall day-to-day variability", terms = c(main_terms, combined_terms, "L2_total")),
  list(group = "Sensitivity", label = "Deaths in first 2 years excluded", terms = main_terms, early = TRUE))
all_terms <- unique(unlist(lapply(specs, `[[`, "terms")))
level2_terms <- c(l2_sd_terms, "L2_total")

#### 2 Model functions ####
# Both methods use svycoxph with the NHANES strata and PSUs; unweighted = all weights 1 (design-based SEs)
designs <- list(
  unweighted = lapply(cohorts, function(x) svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~1, nest = TRUE, data = x$d)),
  weighted = lapply(cohorts, function(x) svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~w_ipw, nest = TRUE, data = x$d)))
cox_formula <- function(terms, event) as.formula(paste0("Surv(permth_exm / 12, ", event, ") ~ ", paste(terms, collapse = " + ")))
cox <- function(d, terms, event, method, ds) {
  svycoxph(cox_formula(terms, event), design = subset(designs[[method]][[ds]], SEQN %in% d$SEQN))
}
# cox.zph is not valid on svycoxph fits, so the PH test uses an ordinary coxph fit (same coefficients as weight 1)
ph_table <- function(d, terms, event) cox.zph(coxph(cox_formula(terms, event), data = d, ties = "efron"))$table
estimate <- function(f, term) {
  b <- unname(coef(f)[term])
  se <- sqrt(unname(vcov(f)[term, term]))
  c(b = b, se = se, p = 2 * pnorm(-abs(b / se)))
}

#### 3 Scores ####
load("./results/clustering/clustering.RData")   # fpca: standard FPCA (6_1)
out_path <- "./results/comparison/"
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
for (ds in names(cohorts)) {
  cohorts[[ds]]$d$tb_FPC4 <- as.numeric(scale(fpca[[ds]]$scores[cohorts[[ds]]$d$SEQN, 4]))
  cohorts[[ds]]$d$FPC2 <- as.numeric(scale(cohorts[[ds]]$d$FPC2))   # hazard ratio per SD
  # Orient standard FPC4 like residual FPC2 (its eigenfunction has the opposite sign), so hazard ratios share a direction
  s4 <- sign(cor(cohorts[[ds]]$d$tb_FPC4, cohorts[[ds]]$d$FPC2, use = "complete.obs"))
  cohorts[[ds]]$d$tb_FPC4 <- s4 * cohorts[[ds]]$d$tb_FPC4
}
designs <- list(
  unweighted = lapply(cohorts, function(x) svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~1, nest = TRUE, data = x$d)),
  weighted = lapply(cohorts, function(x) svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~w_ipw, nest = TRUE, data = x$d)))
demo_terms <- c("Age", "Gender", "Race", "Income", "wear_z")
cosinor_terms <- c("MESOR_z", "Amplitude_z", "cosinor_phase_sin", "cosinor_phase_cos")
# Correlation of standard FPC1-FPC6 and residual FPC1-FPC3 scores with the cosinor parameters; R2 = share explained by all four
cosinor_cor <- bind_rows(lapply(names(cohorts), function(ds) {
  d <- cohorts[[ds]]$d
  S <- cbind(fpca[[ds]]$scores[d$SEQN, 1:6], d[, c("FPC1", "FPC2", "FPC3")])
  colnames(S) <- c(paste0("standard FPC", 1:6), paste0("residual FPC", 1:3))
  bind_rows(lapply(colnames(S), function(v) {
    ok <- complete.cases(S[, v], d[, cosinor_terms])
    data.frame(dataset = ds, score = v, cor(S[ok, v], d[ok, cosinor_terms]),
               R2 = summary(lm(S[ok, v] ~ ., data = d[ok, cosinor_terms]))$r.squared)
  }))
}))
print(cosinor_cor %>% mutate(across(where(is.numeric), ~ round(.x, 2))), row.names = FALSE)

#### 4 Masking: with and without the cosinor parameters ####
masking <- list()
for (ev in c("cvd_death", "mortstat")) for (ds in names(cohorts)) for (m in methods) for (v in c("FPC2", "tb_FPC4")) {
  for (adj in c("demographics", "demographics + cosinor")) {
    d <- cohorts[[ds]]$d
    k <- estimate(cox(d, c(demo_terms, if (adj != "demographics") cosinor_terms, v), ev, m, ds), v)
    masking[[length(masking) + 1]] <- data.frame(outcome = ev, dataset = ds, method = m,
      score = ifelse(v == "FPC2", "residual FPC2", "standard FPC4"), adjustment = adj,
      hr = exp(k[["b"]]), lo = exp(k[["b"]] - 1.96 * k[["se"]]), hi = exp(k[["b"]] + 1.96 * k[["se"]]), p = k[["p"]])
  }
}
masking <- bind_rows(masking)
print(masking %>% mutate(across(c(hr, lo, hi), ~ round(.x, 2)), p = signif(p, 2)), row.names = FALSE)
write.csv(masking, paste0(dir_path, "masking.csv"), row.names = FALSE)

# #### 5 Robustness ####
# swap <- function(terms) c(setdiff(terms, fpc_terms), "tb_FPC4")
# robust_both <- list()
# for (ev in c("cvd_death", "mortstat")) for (ds in names(cohorts)) {
#   d <- cohorts[[ds]]$d
#   d <- d[complete.cases(d[, c("permth_exm", ev, all_terms, "tb_FPC4")]), ]
#   for (m in methods) for (s in specs) {
#     dd <- d
#     if (isTRUE(s$early)) dd <- dd[!(dd$mortstat == 1 & dd$permth_exm < 24), ]
#     for (v in c("FPC2", "tb_FPC4")) {
#       k <- estimate(cox(dd, if (v == "FPC2") s$terms else swap(s$terms), ev, m, ds), v)
#       robust_both[[length(robust_both) + 1]] <- data.frame(outcome = ev, dataset = ds, method = m, set = gsub("\n", " ", s$label),
#         score = ifelse(v == "FPC2", "residual FPC2", "standard FPC4"), n = nrow(dd), deaths = sum(dd[[ev]]),
#         hr = exp(k[["b"]]), lo = exp(k[["b"]] - 1.96 * k[["se"]]), hi = exp(k[["b"]] + 1.96 * k[["se"]]), p = k[["p"]])
#     }
#   }
# }
# robust_both <- bind_rows(robust_both)
# robust_summary <- robust_both %>%
#   group_by(outcome, dataset, method, score) %>%
#   summarise(hr_range = sprintf("%.2f-%.2f", min(hr), max(hr)), significant = paste0(sum(p < 0.05), "/", n()), max_p = max(p), .groups = "drop")
# print(as.data.frame(robust_summary), row.names = FALSE)
