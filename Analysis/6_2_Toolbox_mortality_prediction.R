# Analysis following https://andrew-leroux.github.io/rnhanesdata/articles/vignette_prediction.html
# Mortality prediction with survey-weighted Cox models (svycoxph); deaths from other causes are censored (cause-specific hazard)
#### 0 Packages ####
gc()
library(dplyr)
library(tidyr)
library(survey)
library(survival)
library(ggplot2)
library(patchwork)
library(ggtext)

# cox regression, random survival forest, xgboost
# c-index, AUC, Variable importance, SHAP values
#### 1 Settings ####
outcome <- "cvd_death"
dir_path <- "./results/mortality/"
out_dir <- paste0("./results/mortality_prediction/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
load(paste0(dir_path, "cohorts_weighted.RData"))
# Forward-selection candidates
candidates <- list(Age = "Age", Gender = "Gender", Race = "Race", Income = "Income", WearTime = "wear_z", MESOR = "MESOR_z", Amplitude = "Amplitude_z", Acrophase = c("cosinor_phase_sin", "cosinor_phase_cos"), FPC2 = "FPC2")
accel_candidates <- c("MESOR", "Amplitude", "Acrophase", "FPC2")
base_candidates <- c("Age", "Gender", "Race", "Income", "WearTime")
# Sensitivity analysis: lifestyle and functional variables also compete for selection
lifestyle_candidates <- list(BMI = "BMI", Smoking = "SmokeCigs", Alcohol = "DrinkStatus", Education = "EducationAdult", Mobility = "MobilityProblem")
n_folds <- 10
cohort_label <- c(count = "Hip (2003-2006)", mims = "Wrist (2011-2014)")

#### 2 Helpers ####
fit_cox <- function(design, terms) {
  svycoxph(as.formula(paste("Surv(time, y) ~", paste(terms, collapse = " + "))), design = design)
}
predict_lp <- function(fit, newdata) {
  tt <- delete.response(terms(fit))
  X <- model.matrix(tt, model.frame(tt, newdata, xlev = fit$xlevels))
  as.numeric(X[, names(coef(fit)), drop = FALSE] %*% coef(fit))
}
# Survey-weighted Harrell's C
calc_weighted_C <- function(lp, time, status, weights) {
  unname(concordance(Surv(time, status) ~ lp, weights = weights, reverse = TRUE, timewt = "n")$concordance)
}

# Survey-weighted AUC
calc_weighted_AUC <- function(response, labels, weights) {
  stopifnot(length(response) == length(labels), all(labels %in% c(0, 1)))
  cutpts <- c(Inf, sort(unique(response), decreasing = TRUE))
  f_0 <- weights[labels == 0] / mean(weights[labels == 0])
  f_1 <- weights[labels == 1] / mean(weights[labels == 1])
  x1 <- vapply(cutpts, function(c) mean((response[labels == 1] >= c) * f_1), numeric(1))   # sensitivity
  x2 <- vapply(cutpts, function(c) mean((response[labels == 0] >= c) * f_0), numeric(1))   # 1 - specificity
  n <- length(cutpts)
  sum((x1[-n] + x1[-1]) / 2 * diff(x2))
}

# Cross-validated discrimination:
#   cv_c   = weighted Harrell's C in each test fold, weighted averaged over folds
#   cv_auc = weighted AUC for death within 5 years, from the pooled out-of-fold linear predictors
cv_metrics <- function(d, design, terms) {
  lp <- numeric(nrow(d))
  c_k <- vapply(seq_len(n_folds), function(k) {
    fit <- fit_cox(subset(design, fold != k), terms)
    test <- d$fold == k
    lp[test] <<- predict_lp(fit, d[test, ])
    calc_weighted_C(lp[test], d$time[test], d$y[test], d$w_norm[test])
  }, numeric(1))
  c(cv_c = sum(c_k * tabulate(d$fold, n_folds) / nrow(d)),
    cv_auc = calc_weighted_AUC(lp, as.integer(d$y == 1 & d$time <= 5), d$w_norm))
}

#### 3 Forward selection ####
# Variables are added by CV C-index; selection stops before the first non-increase.
# CV AUC is recorded along the same path. 
forward_select <- function(d, design, cands) {
  inc <- character(0)
  path <- list()
  univariate <- NULL
  while (length(remaining <- setdiff(names(cands), inc)) > 0) {
    step <- bind_rows(lapply(remaining, function(v) {
      terms <- unlist(cands[c(inc, v)], use.names = FALSE)
      fit <- fit_cox(design, terms)
      data.frame(variable = v, t(cv_metrics(d, design, terms)), c_index = calc_weighted_C(fit$linear.predictors, d$time, d$y, d$w_norm))
    }))
    if (length(inc) == 0) univariate <- step %>% arrange(desc(cv_c))   # univariate models
    best <- step[which.max(step$cv_c), ]
    inc <- c(inc, best$variable)
    path[[length(path) + 1]] <- best
  }
  path <- bind_rows(path) %>% mutate(step = row_number(), .before = 1)
  n_sel <- which(diff(path$cv_c) <= 0)[1]
  if (is.na(n_sel)) n_sel <- nrow(path)
  list(path = path, univariate = univariate, selected = path$variable[seq_len(n_sel)])
}

run_selection <- function(candidates) {
forward <- univariate <- final <- final_fit <- discrim <- list()
for (ds in names(cohorts)) {
  d <- cohorts[[ds]]$d
  d$y <- d[[outcome]]
  d$time <- d$permth_exm / 12
  d <- droplevels(d[complete.cases(d[, c("y", "time", unlist(candidates), "w_ipw", "SDMVPSU", "SDMVSTRA")]), ])
  d$w_norm <- d$w_ipw / mean(d$w_ipw)
  # 10 folds stratified by the outcome, so that each test fold has a similar number of deaths
  # set.seed(1244) # andrew-leroux seed
  set.seed(1234) # set the seed so cross-validation results are reproducible
  d$fold <- 0L
  for (ev in 0:1) { idx <- which(d$y == ev); d$fold[idx] <- sample(rep_len(seq_len(n_folds), length(idx))) }
  design <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~w_norm, nest = TRUE, data = d)
  message(ds, ": n = ", nrow(d), ", events (", outcome, ") = ", sum(d$y))

  # Forward selection by CV C-index
  f <- forward_select(d, design, candidates)
  forward[[ds]] <- f$path %>% mutate(dataset = ds, selected = step <= length(f$selected), .before = 1)
  univariate[[ds]] <- f$univariate %>% mutate(dataset = ds, rank = row_number(), .before = 1)

  # Final model (forward selection by CV C-index)
  fit <- fit_cox(design, unlist(candidates[f$selected], use.names = FALSE))
  final_fit[[ds]] <- fit
  final[[ds]] <- data.frame(dataset = ds, term = names(coef(fit)), b = coef(fit), se = sqrt(diag(vcov(fit))), row.names = NULL) %>%
    mutate(p = 2 * pnorm(-abs(b / se)), hr = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se))

  # Discrimination of the base and selected models (same folds)
  sel <- f$path[length(f$selected), ]
  k <- filter(final[[ds]], term == "FPC2")
  fpc2_hr <- if (nrow(k) == 1) sprintf("%.2f (%.2f-%.2f)", k$hr, k$lo, k$hi) else "not selected"
  discrim[[ds]] <- bind_rows(
    data.frame(dataset = ds, model = "Base", variables = paste(base_candidates, collapse = " + "),
               t(cv_metrics(d, design, unlist(candidates[base_candidates], use.names = FALSE))), fpc2_hr = "-"),
    data.frame(dataset = ds, model = "Selected", variables = paste(f$selected, collapse = " + "), cv_c = sel$cv_c, cv_auc = sel$cv_auc, fpc2_hr = fpc2_hr))
}
list(forward = bind_rows(forward), univariate = bind_rows(univariate),
     final = bind_rows(final), final_fit = final_fit, discrim = bind_rows(discrim))
}
main <- run_selection(candidates)
sensitivity <- run_selection(c(candidates, lifestyle_candidates))

#### 4 Plot ####
# CV C-index along the main and sensitivity (+ lifestyle) selection paths
path_cols <- c("Main" = "black", "+ Lifestyle" = "#2360A0")
plot_cohort <- function(ds) {
  paths <- list(`Main` = filter(main$forward, dataset == ds), `+ Lifestyle` = filter(sensitivity$forward, dataset == ds))
  long <- bind_rows(lapply(names(paths), function(nm) {
    p <- paths[[nm]]
    data.frame(step = p$step, line = nm, value = p$cv_c, best = p$step == sum(p$selected))
  })) %>% mutate(line = factor(line, names(path_cols)))
  steps <- sort(unique(long$step))
  var_at <- function(nm, s) {
    v <- paths[[nm]]$variable[paths[[nm]]$step == s]
    if (length(v) == 0) return("")
    if (v %in% accel_candidates) v <- paste0("**", v, "**")
    paste0("<span style='color:", path_cols[[nm]], "'>", v, "</span>")
  }
  labels <- sapply(steps, function(s) paste0("Step ", s, "<br>", var_at("Main", s), "<br>", var_at("+ Lifestyle", s)))
  ggplot(long, aes(step, value, colour = line)) +
    geom_line(linewidth = 0.7) + geom_point(size = 1) +
    geom_point(data = filter(long, best), size = 4, show.legend = FALSE) +
    scale_colour_manual(values = path_cols) +
    scale_x_continuous(breaks = steps, labels = labels) +
    scale_y_continuous(name = "Cross-validated C-index") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_markdown(size = 8, lineheight = 1.2), legend.title = element_blank()) +
    labs(x = NULL, title = cohort_label[[ds]])
}
p_selection <- wrap_plots(lapply(names(cohort_label), plot_cohort), ncol = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
  # plot_annotation(caption = "Forward selection by 10-fold cross-validated, survey-weighted C-index.")
p_selection

#### 5 save ####
print_model <- function(final) print(final %>% transmute(dataset, term, `HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", hr, lo, hi), p = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))), row.names = FALSE)
print_model(main$final)
print_model(sensitivity$final)
# Discrimination table: 10-fold CV, survey-weighted C-index and 5-year AUC
discrimination <- bind_rows(
  main$discrim %>% filter(model == "Base") %>% mutate(model = "Base (demographics)"),
  main$discrim %>% filter(model == "Selected") %>% mutate(model = "Selected (main)"),
  sensitivity$discrim %>% filter(model == "Selected") %>% mutate(model = "Selected (+ lifestyle)")) %>%
  mutate(cohort = cohort_label[dataset], cv_c = round(cv_c, 3), cv_auc = round(cv_auc, 3)) %>%
  arrange(match(dataset, names(cohort_label))) %>%
  select(cohort, model, variables, `CV C-index` = cv_c, `CV 5-year AUC` = cv_auc, `FPC2 HR (95% CI)` = fpc2_hr)
print(discrimination, row.names = FALSE)
write.csv(discrimination, paste0(out_dir, "discrimination.csv"), row.names = FALSE)
saveRDS(main$final_fit, paste0(out_dir, "final_model_fit.rds"))
