# Analysis following https://andrew-leroux.github.io/rnhanesdata/articles/vignette_prediction.html
# Mortality prediction with survey-weighted logistic models (svyglm, quasibinomial)
#### 0 Packages ####
gc()
library(dplyr)
library(tidyr)
library(survey)
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
candidates <- list(Age = "Age", Gender = "Gender", Race = "Race", Income = "Income", MESOR = "MESOR_z", Amplitude = "Amplitude_z", Acrophase = c("cosinor_phase_sin", "cosinor_phase_cos"), FPC2 = "FPC2")
accel_candidates <- c("MESOR", "Amplitude", "Acrophase", "FPC2")
# Sensitivity analysis: lifestyle and functional variables also compete for selection
lifestyle_candidates <- list(BMI = "BMI", Smoking = "SmokeCigs", Alcohol = "DrinkStatus", Education = "EducationAdult", Mobility = "MobilityProblem")
n_folds <- 10
# fold_seeds <- c(1244) # andrew-leroux seed
fold_seeds <- c(1234)
cohort_label <- c(count = "Hip (2003-2006)", mims = "Wrist (2011-2014)")

#### 2 Helpers ####
fit_glm <- function(design, terms) {
  svyglm(as.formula(paste("y ~", paste(terms, collapse = " + "))), design = design, family = quasibinomial())
}
# Complex survey AIC (Lumley and Scott 2015); 
svy_aic <- function(fit, k = 2) unname(extractAIC(fit, k = k)[2])
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

# Cross-validated weighted AUC
cv_auc <- function(d, design, terms) {
  auc_k <- vapply(seq_len(n_folds), function(k) {
    fit <- fit_glm(subset(design, fold != k), terms)
    test <- d[d$fold == k, ]
    lp <- as.numeric(predict(fit, newdata = test, type = "link", se.fit = FALSE))
    calc_weighted_AUC(lp, test$y, test$w_norm)
  }, numeric(1))
  sum(auc_k * tabulate(d$fold, n_folds) / nrow(d))
}

#### 3 Forward selection ####
forward_select <- function(d, design, cands, criterion = "cv_auc") {
  better <- if (criterion == "cv_auc") which.max else which.min
  inc <- character(0)
  path <- list()
  univariate <- NULL
  while (length(remaining <- setdiff(names(cands), inc)) > 0) {
    step <- bind_rows(lapply(remaining, function(v) {
      terms <- unlist(cands[c(inc, v)], use.names = FALSE)
      fit <- fit_glm(design, terms)
      data.frame(variable = v, cv_auc = cv_auc(d, design, terms), auc = calc_weighted_AUC(fit$linear.predictors, d$y, d$w_norm), AIC = svy_aic(fit), EPIC = svy_aic(fit, k = 4))
    }))
    if (all(is.na(step$AIC))) break
    if (length(inc) == 0) univariate <- step %>% arrange(desc(cv_auc))   # univariate models
    best <- step[better(step[[criterion]]), ]
    inc <- c(inc, best$variable)
    path[[length(path) + 1]] <- best
  }
  path <- bind_rows(path) %>% mutate(criterion = criterion, step = row_number(), .before = 1)
  change <- diff(path[[criterion]])
  n_sel <- which(if (criterion == "cv_auc") change < 0 else change > 0)[1]
  if (is.na(n_sel)) n_sel <- nrow(path)
  list(path = path, univariate = univariate, selected = path$variable[seq_len(n_sel)])
}

run_selection <- function(candidates) {
forward <- forward_ic <- univariate <- final <- final_fit <- list()
for (ds in names(cohorts)) {
  d <- cohorts[[ds]]$d
  d$y <- d[[outcome]]
  d <- droplevels(d[complete.cases(d[, c("y", unlist(candidates), "w_ipw", "SDMVPSU", "SDMVSTRA")]), ])
  d$w_norm <- d$w_ipw / mean(d$w_ipw)
  # 10 folds stratified by the outcome, so that each test fold has a similar number of deaths
  # set.seed(1244) # andrew-leroux seed
  set.seed(1234) # set the seed so cross-validation results are reproducible
  d$fold <- 0L
  for (ev in 0:1) { idx <- which(d$y == ev); d$fold[idx] <- sample(rep_len(seq_len(n_folds), length(idx))) }
  design <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~w_norm, nest = TRUE, data = d)
  message(ds, ": n = ", nrow(d), ", events (", outcome, ") = ", sum(d$y))

  # Forward selection by CV AUC (main), EPIC and AIC
  f <- forward_select(d, design, candidates)
  forward[[ds]] <- f$path %>% mutate(dataset = ds, selected = step <= length(f$selected), .before = 1)
  univariate[[ds]] <- f$univariate %>% mutate(dataset = ds, rank = row_number(), .before = 1)
  for (cr in c("EPIC", "AIC")) {
    fc <- forward_select(d, design, candidates, cr)
    forward_ic[[paste(ds, cr)]] <- fc$path %>% mutate(dataset = ds, selected = step <= length(fc$selected), .before = 1)
  }

  # Final model (forward selection by CV AUC)
  fit <- fit_glm(design, unlist(candidates[f$selected], use.names = FALSE))
  final_fit[[ds]] <- fit
  cf <- summary(fit, df.resid = Inf)$coefficients
  final[[ds]] <- data.frame(dataset = ds, term = rownames(cf), b = cf[, 1], se = cf[, 2], p = cf[, 4], row.names = NULL) %>%
    mutate(or = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se))
}
list(forward = bind_rows(forward), forward_ic = bind_rows(forward_ic), univariate = bind_rows(univariate),
     final = bind_rows(final), final_fit = final_fit)
}
main <- run_selection(candidates)
sensitivity <- run_selection(c(candidates, lifestyle_candidates))

#### 4 Plot ####
path_cols <- c("CV AUC" = "black", "EPIC" = "#2360A0", "AIC" = "#3B8B3B")
plot_criteria <- function(ds, res) {
  paths <- list(cv_auc = filter(res$forward, dataset == ds), EPIC = filter(res$forward_ic, dataset == ds, criterion == "EPIC"), AIC = filter(res$forward_ic, dataset == ds, criterion == "AIC"))
  paths <- lapply(paths, function(p) filter(p, !is.na(AIC)))
  r_ic <- range(c(paths$EPIC$EPIC, paths$AIC$AIC))
  r_auc <- range(paths$cv_auc$cv_auc)
  b <- diff(r_ic) / diff(r_auc)
  a <- r_ic[1] - b * r_auc[1]
  long <- bind_rows(lapply(names(paths), function(cr) {
    p <- paths[[cr]]
    data.frame(step = p$step, variable = p$variable, line = names(path_cols)[match(cr, names(paths))], value = if (cr == "cv_auc") a + b * p$cv_auc else p[[cr]], best = p$step == sum(p$selected))
  })) %>% mutate(line = factor(line, names(path_cols)))
 
  steps <- sort(unique(long$step))
  var_at <- function(cr, s) {
    v <- paths[[cr]]$variable[paths[[cr]]$step == s]
    if (length(v) == 0) return("")
    if (v %in% accel_candidates) v <- paste0("**", v, "**")
    paste0("<span style='color:", path_cols[match(cr, names(paths))], "'>", v, "</span>")
  }
  labels <- sapply(steps, function(s) paste0("Step ", s, "<br>", var_at("cv_auc", s), "<br>", var_at("EPIC", s), "<br>", var_at("AIC", s)))
  ggplot(long, aes(step, value, colour = line)) +
    geom_line(linewidth = 0.7) + geom_point(size = 1) +
    geom_point(data = filter(long, best), size = 4, show.legend = FALSE) +
    scale_colour_manual(values = path_cols) +
    scale_x_continuous(breaks = steps, labels = labels) +
    scale_y_continuous(name = "Survey AIC / EPIC", sec.axis = sec_axis(~ (. - a) / b, name = "Cross-validated AUC")) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_markdown(size = 8, lineheight = 1.2), legend.title = element_blank()) +
    labs(x = NULL, title = cohort_label[[ds]])
}

plot_selection <- function(res) {
  wrap_plots(lapply(names(cohort_label), plot_criteria, res = res), ncol = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
    # plot_annotation(caption = paste("Forward selection by each criterion (10-fold cross-validated, survey-weighted AUC; survey AIC; EPIC = AIC with penalty 4).", sep = "\n"))  
}
p_main <- plot_selection(main)
p_sensitivity <- plot_selection(sensitivity)
p_main
p_sensitivity

#### 5 save ####
print_model <- function(final) print(final %>% transmute(dataset, term, `OR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", or, lo, hi), p = ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))), row.names = FALSE)
print_model(main$final)
print_model(sensitivity$final)
saveRDS(main$final_fit, paste0(out_dir, "final_model_fit.rds"))
