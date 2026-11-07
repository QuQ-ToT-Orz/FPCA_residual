#### 0 Packages ####
gc()
library(dplyr)
library(tidyr)
library(survival)
library(survey)
library(haven)
library(splines)
library(rnhanesdata)

#### 1 Survey weights ####
# rnhanesdata::reweight_accel()
# Survey weights for an accelerometry analytic sample that is a subset of the NHANES MEC-examined sample.
# MEC examination weights divided by the estimated probability of inclusion in the analytic sample (inverse probability of selection), truncated at the 99th percentile.
expanded <- TRUE   # FALSE for the default run
out_suffix <- if (expanded) "_expanded" else ""
dir_path <- paste0("./results/mortality", out_suffix, "/")
load(paste0(dir_path, "cohorts.RData"))
truncate <- 0.99
with_missing <- function(x) { x <- as.character(x); x[is.na(x)] <- "Missing"; factor(x) }

# MEC-examined participants of the same cycles
load_mec_covariates <- function(dataset) {
  e <- new.env()
  if (dataset == "mims") {
    load("./data/mims/Covariate_G.RData", e)
    load("./data/mims/Covariate_H.RData", e)
    covariates <- bind_rows(e$Covariate_G, e$Covariate_H)
  } else {
    data(list = c("Covariate_C", "Covariate_D"), package = "rnhanesdata", envir = e)
    covariates <- bind_rows(e$Covariate_C, e$Covariate_D)
  }
  income <- bind_rows(lapply(if (dataset == "mims") c("DEMO_G", "DEMO_H") else c("DEMO_C", "DEMO_D"), function(f) {
    read_xpt(paste0("./data/nhanes_demo/", f, ".xpt"))[, c("SEQN", "INDFMPIR")]
  })) %>%
    transmute(SEQN = as.character(as.integer(SEQN)),
              income_cat = cut(INDFMPIR, c(-Inf, 1.3, 3.5, Inf), right = FALSE, labels = c("Low", "Middle", "High")))
  # Covariate files keep the original race levels; merge the Hispanic groups as in 1_data_*
  covariates %>%
    mutate(SEQN = as.character(SEQN), Race = as.character(Race),
           Race = ifelse(Race %in% c("Mexican American", "Other Hispanic"), "Hispanic", Race)) %>%
    left_join(income, by = "SEQN")
}

#### 2 Selection model and weights ####
# The model is fitted with the MEC weights and the NHANES design (survey::svyglm), so that the inclusion probabilities refer to the examined population.
weights <- list()
for (ds in names(cohorts)) {
  d <- cohorts[[ds]]$d
  min_age <- floor(min(d$Age)) # 20
  max_age <- floor(max(d$Age))
  population <- load_mec_covariates(ds) %>%
    filter(WTMEC2YR > 0, RIDAGEYR >= min_age, RIDAGEYR <= max_age) %>%
    mutate(included = as.integer(SEQN %in% d$SEQN), cycle = factor(SDDSRVYR),
           across(c(Gender, Race, income_cat, EducationAdult, BMI_cat, SmokeCigs, DrinkStatus, Diabetes, CHD, CHF, Stroke, Cancer, MobilityProblem), with_missing))
  design <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC2YR, nest = TRUE, data = population)
  fit <- svyglm(included ~ cycle + ns(RIDAGEYR, 4) + Gender + Race + income_cat + EducationAdult + BMI_cat + SmokeCigs + DrinkStatus + Diabetes + CHD + CHF + Stroke + Cancer + MobilityProblem, design = design, family = quasibinomial())
  population$p <- as.numeric(fitted(fit))
  p <- population$p[match(d$SEQN, population$SEQN)]
  w <- d$WTMEC2YR / 2 / p
  cap <- quantile(w, truncate)
  weights[[ds]] <- data.frame(dataset = ds, SEQN = d$SEQN, p_included = p, w_ipw = pmin(w, cap))
}
weights <- bind_rows(weights)
for (ds in names(cohorts)) {
  cohorts[[ds]]$d <- cohorts[[ds]]$d %>%
    left_join(weights %>% filter(dataset == ds) %>% select(SEQN, w_ipw), by = "SEQN")
}

save(cohorts, file = paste0(dir_path, "cohorts_weighted.RData"))

#### 3 Survival models ####
# Model 1: age, sex, race/ethnicity, income, cosinor MESOR, amplitude and acrophase (weekly fit), FPC1-FPC3, for all-cause, cardiovascular and cancer mortality (cause-specific: other deaths censored).
# Testing hierarchy: the hip cohort (discovery) searches FPC1-FPC3 for all-cause mortality, with Benjamini-Hochberg q-values; the wrist cohort (validation) and cause-specific mortality test the component carried forward (FPC2) as a single test. FPC2 is then pooled across cohorts (descriptive).
# Methods: svycoxph with NHANES strata and PSUs, unweighted (all weights 1) and survey-weighted (IPW weights).
base_terms <- c("Age", "Gender", "Race", "Income", "MESOR_z", "Amplitude_z", "cosinor_phase_sin", "cosinor_phase_cos")
fpc_terms <- c("FPC1", "FPC2", "FPC3")
outcomes <- c(mortstat = "All-cause", cvd_death = "Cardiovascular", cancer_death = "Cancer")
methods <- c("unweighted", "weighted")

# Both methods use svycoxph with the NHANES strata and PSUs; unweighted = all weights 1 (design-based SEs)
designs <- list(
  unweighted = lapply(cohorts, function(x) svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~1, nest = TRUE, data = x$d)),
  weighted = lapply(cohorts, function(x) svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~w_ipw, nest = TRUE, data = x$d)))
cox_formula <- function(terms, event) as.formula(paste0("Surv(permth_exm / 12, ", event, ") ~ ", paste(terms, collapse = " + ")))
cox <- function(d, terms, event, method, ds) {
  svycoxph(cox_formula(terms, event), design = subset(designs[[method]][[ds]], SEQN %in% d$SEQN))
}
# cox.zph is not valid on svycoxph fits, so the PH test uses an ordinary coxph fit (same coefficients as weight 1)
ph_table <- function(d, terms, event) cox.zph(coxph(cox_formula(terms, event), data = d))$table
estimate <- function(f, term) {
  b <- unname(coef(f)[term])
  se <- sqrt(unname(vcov(f)[term, term]))
  c(b = b, se = se, p = 2 * pnorm(-abs(b / se)))
}

#### 4 Model 1 in each cohort ####
terms <- c(base_terms, fpc_terms)
results <- list()
for (ev in names(outcomes)) for (ds in names(cohorts)) for (m in methods) {
  d <- cohorts[[ds]]$d
  d <- d[complete.cases(d[, c("permth_exm", ev, terms)]), ]
  f <- cox(d, terms, ev, m, ds)
  ph <- if (m == "unweighted") ph_table(d, terms, ev) else NULL   # proportional hazards, unweighted fit
  for (tm in c(fpc_terms, "MESOR_z", "Amplitude_z")) {
    k <- estimate(f, tm)
    results[[length(results) + 1]] <- data.frame(
      outcome = ev, dataset = ds, method = m, term = tm, n = nrow(d), deaths = sum(d[[ev]]),
      b = k[["b"]], se = k[["se"]], p = k[["p"]],
      ph_p = if (is.null(ph)) NA_real_ else ph[tm, "p"], # proportional-hazards diagnostics
      ph_global_p = if (is.null(ph)) NA_real_ else ph["GLOBAL", "p"])
  }
}
results <- bind_rows(results) %>%
  mutate(role = ifelse(dataset == "count", "discovery (hip)", "validation (wrist)"),
         hr = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se)) %>%
  group_by(outcome, dataset, method) %>%
  mutate(q = { q <- rep(NA_real_, n())
               if (dataset[1] == "count" && outcome[1] == "mortstat") q[term %in% fpc_terms] <- p.adjust(p[term %in% fpc_terms], "BH")
               q }) %>%
  ungroup()

#### 5 Pooled FPC2 across cohorts ####
# Fixed-effect inverse-variance pooling; heterogeneity: Cochran's Q.
pooled <- results %>%
  filter(term == "FPC2") %>%
  group_by(outcome, method, term) %>%
  summarise(het_p = pchisq(sum((b - sum(b / se^2) / sum(1 / se^2))^2 / se^2), 1, lower.tail = FALSE), # Heterogeneity
            b = sum(b / se^2) / sum(1 / se^2), se = sqrt(1 / sum(1 / se^2)),
            n = sum(n), deaths = sum(deaths), .groups = "drop") %>%
  mutate(dataset = "pooled", role = "pooled", p = 2 * pnorm(-abs(b / se)),
         hr = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se))
results <- bind_rows(results, pooled)

#### 6 Save ####
pfmt <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
ph_unweighted <- results %>%
  filter(method == "unweighted", dataset != "pooled") %>%
  select(outcome, dataset, term, ph_test = ph_p)
res_table <- function(m) {
  results %>%
    filter(method == m, term %in% fpc_terms) %>%
    select(-ph_p) %>%
    left_join(ph_unweighted, by = c("outcome", "dataset", "term")) %>%
    mutate(hr_ci = sprintf("%.2f (%.2f-%.2f)", hr, lo, hi),
           p_txt = ifelse(is.na(q), pfmt(p), paste0(pfmt(p), " (q=", pfmt(q), ")")),
           ph_txt = ifelse(is.na(ph_test), "", pfmt(ph_test)),
           het_txt = ifelse(is.na(het_p), "", pfmt(het_p))) %>%
    select(outcome, term, dataset, hr_ci, p_txt, ph_txt, het_txt) %>%
    pivot_wider(names_from = dataset, values_from = c(hr_ci, p_txt, ph_txt, het_txt), values_fill = "") %>%
    arrange(factor(outcome, names(outcomes)), term) %>%
    transmute(Outcome = ifelse(duplicated(outcome), "", outcomes[outcome]), Component = term,
    `Hip HR (95% CI)` = hr_ci_count, `Hip p` = p_txt_count, `Hip PH p` = ph_txt_count,
    `Wrist HR (95% CI)` = hr_ci_mims, `Wrist p` = p_txt_mims, `Wrist PH p` = ph_txt_mims,
    `Pooled HR (95% CI)` = hr_ci_pooled, `Pooled p` = p_txt_pooled, `Heterogeneity p` = het_txt_pooled)
}
# Participants and deaths in each model
sample_size <- results %>%
  filter(method == "weighted", term == "FPC2", dataset != "pooled") %>%
  transmute(Outcome = outcomes[outcome], dataset, n_deaths = sprintf("%d / %d", deaths, n)) %>%
  pivot_wider(names_from = dataset, values_from = n_deaths) %>%
  rename(`Hip deaths / n` = count, `Wrist deaths / n` = mims)
table_weighted <- res_table("weighted")       # main: survey-weighted (IPW)
table_unweighted <- res_table("unweighted")   # supplement
print(as.data.frame(table_unweighted), row.names = FALSE)
print(as.data.frame(table_weighted), row.names = FALSE)
print(as.data.frame(sample_size), row.names = FALSE)
#         Outcome Hip deaths / n Wrist deaths / n
#       All-cause     463 / 2615       312 / 3208
#  Cardiovascular     112 / 2615       106 / 3208
#          Cancer     101 / 2615        77 / 3208
# expanded
#         Outcome Hip deaths / n Wrist deaths / n
#       All-cause     507 / 2917       362 / 3605
#  Cardiovascular     120 / 2917       125 / 3605
#          Cancer     109 / 2917        91 / 3605
write.csv(results, paste0(dir_path, "fpc_mortality.csv"), row.names = FALSE)
