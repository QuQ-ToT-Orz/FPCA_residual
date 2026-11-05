# Apply the fitted FACE score to new participants, 
# using training means, eigenfunctions,and shrinkage
# No new model is fit.
#### 0 Packages ####
gc()
library(dplyr)
library(ggplot2)
library(ActCR)
library(survey)

cosinor_fit <- "weekly"
dir_path <- "./results/score_reference/"
dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
time <- 1:1440

#### 1 Reference objects from the fitted cohorts ####
# projection: P_ik = Σ_t (Y_i(t) − μ(t)) φ_k(t) = ξ_ik + e_ik,    e_ik ~ (0, σ²)
# shrinkage c_k: ξ̂_ik = E[ξ_ik | P_ik] = λ_k / (λ_k + σ²) · P_ik # shrinkage <- f$evalues[2] /(f$evalues[2] + f$sigma2)
# score = shrinkage x projection # score <- sum((residual - f$mu) * f$efunctions[, 2]) * shrinkage
# c_k is one constant per component, shared by everyone
# λ_k are the eigenvalues, σ² is minute-level noise
load("./results/mortality/cohorts_weighted.RData") 
build_reference <- function(ds) {
  suffix <- if (ds == "mims") "new" else "old"
  e <- new.env()
  load(paste0("./results/cosinor_residual_fpca/", ds, "/residual_fpca_", cosinor_fit, "_", suffix, ".RData"), envir = e)
  f <- e$fpca_fit
  phi <- f$efunctions[, 2]
  # Manual calculation
  # projection <- drop(sweep(e$subject_residual[rownames(f$scores), ], 2, f$mu) %*% phi)
  # shrinkage <- sum(f$scores[, 2] * projection) / sum(projection^2)
  shrinkage <- f$shrinkage[2]

  sign <- sign(sum(phi * cohorts[[ds]]$efunctions[, 2]))
  scores <- sign * f$scores[cohorts[[ds]]$d$SEQN, 2]
  list(dataset = ds, unit = if (ds == "mims") "MIMS units (wrist)" else "activity counts (hip)",
       mu = f$mu, efunction = sign * phi, shrinkage = shrinkage,
       score_mean = mean(scores), score_sd = sd(scores), cohort_scores = scores)
}
reference <- lapply(c(count = "count", mims = "mims"), build_reference)
saveRDS(reference, paste0(dir_path, "score_reference.rds"))

#### 2 Scoring function ####
score_fpc <- function(minute_data, dataset = c("count", "mims"), reference_path = paste0(dir_path, "score_reference.rds")) {
  dataset <- match.arg(dataset)
  ref <- readRDS(reference_path)[[dataset]]
  x <- as.matrix(minute_data)
  x[is.na(x) | x < 0] <- 0
  x <- asinh(x)
  p <- as.data.frame(ActCosinor(c(t(x)), window = 1, export_ts = FALSE)$params)
  acrotime <- p$acrotime * 60
  cosinor <- p$mes + p$amp * cos(2 * pi * (time - acrotime) / 1440)
  residual <- colMeans(sweep(x, 2, cosinor))
  score <- sum(ref$efunction * (residual - ref$mu)) * ref$shrinkage
  list(
    cosinor = data.frame(MESOR = p$mes, amplitude = p$amp, acrophase_minute = acrotime, acrophase = sprintf("%02d:%02d", floor(acrotime / 60) %% 24, round(acrotime %% 60) %% 60)),
    FPC2 = data.frame(score = score, cohort_sd_units = (score - ref$score_mean) / ref$score_sd, cohort_percentile = 100 * mean(ref$cohort_scores <= score)),
    n_days = nrow(x),
    curves = data.frame(minute = time, residual = residual,
                        fpc2_part = ref$mu + score * ref$efunction, cohort_mean = ref$mu))
}

# Residual curve, the cohort mean residual plus the participant's FPC2 component, and the cohort mean residual
plot_fpc <- function(result) {
  d <- tidyr::pivot_longer(result$curves, -minute, names_to = "curve")
  d$curve <- factor(d$curve, c("residual", "fpc2_part", "cohort_mean"), c("Participant residual", "Cohort mean + FPC2 component", "Cohort mean residual"))
  ggplot(d, aes(minute, value, colour = curve, linewidth = curve)) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_line() +
    scale_colour_manual(values = c("grey60", "#B74F27", "#2360A0")) +
    scale_linewidth_manual(values = c(0.3, 0.9, 0.7)) +
    scale_x_continuous(breaks = c(1, 241, 481, 721, 961, 1201, 1440), labels = c("12 AM", "4 AM", "8 AM", "12 PM", "4 PM", "8 PM", "12 AM")) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom") +
    labs(x = NULL, y = "Cosinor residual (asinh units)", colour = NULL, linewidth = NULL,
         subtitle = with(result$FPC2, sprintf("FPC2 = %.2f SD (%.0fth percentile)", cohort_sd_units, cohort_percentile)))
}

# Two example participants
example_ids <- c(count = "30858", mims = "69831")
examples <- check <- list()
for (ds in names(example_ids)) {
    i <- example_ids[[ds]]
    suffix <- if (ds == "mims") "new" else "old"
    load(paste0("./data/", ds, "/Act_Analysis_", suffix, ".RData"))
    Act_Analysis <- Act_Analysis %>% arrange(SEQN, WEEKDAY)
    r <- score_fpc(Act_Analysis[as.character(Act_Analysis$SEQN) == i, paste0("MIN", 1:1440)], ds)
    rm(Act_Analysis)
    examples[[paste(ds, i)]] <- r
    check[[ds]] <- data.frame(dataset = ds, SEQN = i, rescored = r$FPC2$cohort_sd_units, saved = cohorts[[ds]]$d$FPC2[cohorts[[ds]]$d$SEQN == i])
    print(cbind(dataset = ds, SEQN = i, r$cosinor, r$FPC2), row.names = FALSE)
    print(plot_fpc(r) + ggtitle(paste0(if (ds == "mims") "Wrist MIMS" else "Hip counts", ", SEQN ", i)))
}
print(bind_rows(check) %>% mutate(across(c(rescored, saved), ~ round(.x, 4))), row.names = FALSE)
# dataset  SEQN    MESOR amplitude acrophase_minute acrophase    score cohort_sd_units cohort_percentile
#    count 30858 1.223872  1.591262         932.8335     15:33 15.49183        1.546237          93.00191
#  dataset  SEQN   MESOR amplitude acrophase_minute acrophase   score cohort_sd_units cohort_percentile
#     mims 69831 1.79084  1.776485          843.367     14:03 1.66498       0.2586919          63.31047

#### 3 Mortality prediction based on forward selection ####
final_fit <- readRDS("./results/mortality_prediction/final_model_fit.rds")
predict_risk <- function(d, ds) plogis(as.numeric(predict(final_fit[[ds]], newdata = d, type = "link", se.fit = FALSE)))
person_vars <- function(dataset) setdiff(all.vars(formula(final_fit[[dataset]]))[-1],
  c("MESOR_z", "Amplitude_z", "cosinor_phase_sin", "cosinor_phase_cos", "FPC2"))
predict_mortality <- function(result, dataset, ...) {
  cohort <- cohorts[[dataset]]$d
  person <- list(...)
  needed <- person_vars(dataset)
  if (!setequal(names(person), needed)) stop("Give exactly the variables in the ", dataset, " model: ", paste(needed, collapse = ", "))
  for (v in Filter(function(v) is.factor(cohort[[v]]), needed)) stopifnot(person[[v]] %in% levels(cohort[[v]]))
  new <- data.frame(
    person,
    MESOR_z = (result$cosinor$MESOR - mean(cohort$cosinor_mes)) / sd(cohort$cosinor_mes),
    Amplitude_z = (result$cosinor$amplitude - mean(cohort$cosinor_amp)) / sd(cohort$cosinor_amp),
    cosinor_phase_sin = sin(2 * pi * result$cosinor$acrophase_minute / 1440),
    cosinor_phase_cos = cos(2 * pi * result$cosinor$acrophase_minute / 1440),
    FPC2 = result$FPC2$cohort_sd_units)

  if ("Age" %in% needed && (person$Age < min(cohort$Age) || person$Age > max(cohort$Age))) warning("Age is outside the range of the ", dataset, " cohort.")
  p <- predict(final_fit[[dataset]], newdata = new, type = "link", se.fit = TRUE)
  lp <- as.numeric(p)
  se <- as.numeric(SE(p))
  risk <- plogis(lp)
  cohort_risk <- predict_risk(cohort, dataset)
  data.frame(risk = risk, lo = plogis(lp - 1.96 * se), hi = plogis(lp + 1.96 * se), 
  cohort_mean_risk = weighted.mean(cohort_risk, cohort$w_ipw), 
  cohort_percentile = 100 * mean(cohort_risk <= risk), row.names = NULL)
}

# Two example participants
risk_examples <- bind_rows(lapply(names(example_ids), function(ds) {
  i <- example_ids[[ds]]
  x <- cohorts[[ds]]$d[cohorts[[ds]]$d$SEQN == i, ]
  r <- do.call(predict_mortality, c(list(examples[[paste(ds, i)]], ds), as.list(x[person_vars(ds)])))
  data.frame(dataset = ds, SEQN = i, r, risk_saved = predict_risk(x, ds))
}))
print(risk_examples %>% mutate(across(where(is.double), ~ round(.x, 4))), row.names = FALSE)
