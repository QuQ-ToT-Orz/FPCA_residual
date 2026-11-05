# Builds the app's reference objects and example files from the results of Analysis/ scripts 1-6.
# Run from the project root: Rscript Shiny_tool/build_reference.R
library(dplyr)
library(survey)

load("./results/mortality/cohorts_weighted.RData")   # cohorts
load("./results/clustering/clustering.RData")        # fpca, clustering, cluster_summary, profile_curves
score_reference <- readRDS("./results/score_reference/score_reference.rds")
final_fit <- readRDS("./results/mortality_prediction/final_model_fit.rds")
cox <- read.csv("./results/mortality/fpc_mortality.csv") %>%
  filter(term == "FPC2", method == "weighted", outcome == "cvd_death")
labels <- c(count = "Hip counts, NHANES 2003-2006", mims = "Wrist MIMS, NHANES 2011-2014")
example_name <- c(count = "hip", mims = "wrist")
suffix <- c(count = "old", mims = "new")

# Keeps what predict.svyglm needs and drops participant-level data, so the app can be shared
slim_fit <- function(f) {
  f[c("data", "model", "y", "x", "residuals", "fitted.values", "linear.predictors", "weights", "prior.weights", "effects", "qr", "R", "survey.design")] <- NULL
  attr(f$terms, ".Environment") <- globalenv()
  environment(f$formula) <- globalenv()
  f
}

build_reference <- function(ds) {
  d <- cohorts[[ds]]$d
  d$cluster <- clustering[[ds]]$cluster[d$SEQN]
  wmean <- function(x, w = d$w_ipw) weighted.mean(x, w, na.rm = TRUE)
  summarise_group <- function(g) summarise(g,
    Participants = n(), `Survey-weighted share, %` = 100 * sum(w_ipw) / sum(d$w_ipw),
    `Age, years` = weighted.mean(Age, w_ipw), `Female, %` = 100 * weighted.mean(Gender == "Female", w_ipw),
    `BMI, kg/m2` = weighted.mean(BMI, w_ipw, na.rm = TRUE),
    MESOR = sprintf("%.2f", weighted.mean(cosinor_mes, w_ipw)),
    Amplitude = sprintf("%.2f", weighted.mean(cosinor_amp, w_ipw)),
    # Circular mean of the acrophase, as clock time
    Acrophase = {
      minute <- (atan2(weighted.mean(cosinor_phase_sin, w_ipw), weighted.mean(cosinor_phase_cos, w_ipw)) * 1440 / (2 * pi)) %% 1440
      sprintf("%02d:%02d", floor(minute / 60), floor(minute %% 60))
    },
    `Mean FPC2, SD` = sprintf("%.2f", weighted.mean(FPC2, w_ipw)))
  profile <- bind_rows(
    d %>% group_by(Cluster = as.character(cluster)) %>% summarise_group(),
    d %>% mutate(Cluster = "All") %>% group_by(Cluster) %>% summarise_group())
  fit <- final_fit[[ds]]
  risk <- plogis(as.numeric(predict(fit, newdata = d, type = "link", se.fit = FALSE)))
  terms <- setdiff(names(coef(fit)), "(Intercept)")
  terms <- sub("^(Gender|Income).*", "\\1", sub("cosinor_phase_.*", "Acrophase", sub("_z$", "", terms)))
  r <- score_reference[[ds]]
  list(
    label = labels[[ds]],
    followup = median(d$permth_exm / 12),
    levels = list(gender = levels(d$Gender), income = levels(d$Income)),
    cosinor = c(mes_mean = mean(d$cosinor_mes), mes_sd = sd(d$cosinor_mes),
                amp_mean = mean(d$cosinor_amp), amp_sd = sd(d$cosinor_amp)),
    fpc2 = list(mu = r$mu, efunction = r$efunction, shrinkage = r$shrinkage,
                mean = r$score_mean, sd = r$score_sd, scores = sort(unname(r$cohort_scores)),
                hr = unlist(cox[cox$dataset == ds, c("hr", "lo", "hi")])),
    cluster = list(mu = fpca[[ds]]$mu, efunctions = fpca[[ds]]$efunctions, shrinkage = fpca[[ds]]$shrinkage,
                   centers = clustering[[ds]]$centers, profile = as.data.frame(profile),
                   silhouette = cluster_summary$silhouette[cluster_summary$dataset == ds],
                   curves = filter(profile_curves, dataset == ds, scale == "asinh", curve == "Observed")),
    mortality = list(fit = slim_fit(fit), risks = sort(risk), mean = wmean(risk), terms = unique(terms)))
}
reference <- lapply(c(count = "count", mims = "mims"), build_reference)
saveRDS(reference, "./Shiny_tool/reference.rds")

# Examples: the participant closest to each cluster centre
dir.create("./Shiny_tool/examples", showWarnings = FALSE)
examples <- bind_rows(lapply(names(reference), function(ds) {
  scores <- fpca[[ds]]$scores
  cl <- clustering[[ds]]$cluster
  distance <- sqrt(rowSums((scores - clustering[[ds]]$centers[cl, ])^2))
  ids <- sapply(sort(unique(cl)), function(k) names(which.min(ifelse(cl == k, distance, Inf))))
  e <- new.env()
  load(paste0("./data/", ds, "/Act_Analysis_", suffix[ds], ".RData"), e)
  columns <- c("SEQN", "WEEKDAY", paste0("MIN", 1:1440))
  d <- cohorts[[ds]]$d[match(ids, cohorts[[ds]]$d$SEQN), ]
  file <- paste0(example_name[ds], "-cluster-", cl[ids])
  for (i in seq_along(ids)) write.csv(e$Act_Analysis[e$Act_Analysis$SEQN == ids[i], columns],
                                       paste0("./Shiny_tool/examples/", file[i], ".csv"), row.names = FALSE)
  data.frame(file = paste0(file, ".csv"), cohort = ds, SEQN = ids, cluster = cl[ids], age = d$Age, gender = d$Gender, race = d$Race, income = d$Income)
}))
write.csv(examples, "./Shiny_tool/examples/examples.csv", row.names = FALSE)
