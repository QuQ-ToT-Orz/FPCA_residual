#### 0 Packages ####
gc()
library(dplyr)
library(tidyr)
library(ggplot2)

#### 1 Load the data ####
dataset <- "mims"
dataset <- "count"
cosinor_fit <- "weekly"   # weekly (default) or "daily"
expanded <- FALSE   # FALSE for the default run
out_suffix <- if (expanded) "_expanded" else ""
suffix <- paste0(if (dataset == "mims") "new" else "old", out_suffix)
unit <- if (dataset == "mims") "asinh(MIMS)" else "asinh(activity count)"
dir_path <- paste0("./results/cosinor_residual_fpca/", dataset, "/")
load(paste0(dir_path, "residual_fpca_", cosinor_fit, "_", suffix, ".RData"))
load(paste0("./data/", dataset, "/Act_Analysis_asinh_", suffix, ".RData"))
Act_Analysis_asinh <- Act_Analysis_asinh %>% arrange(SEQN, WEEKDAY)
observed <- as.matrix(Act_Analysis_asinh[, paste0("MIN", 1:1440)])
day_id <- as.character(Act_Analysis_asinh$SEQN)
subject_observed <- t(sapply(split(seq_len(nrow(observed)), day_id),
                             function(inx) colMeans(observed[inx, , drop = FALSE])))
rm(observed, Act_Analysis_asinh)

time <- 1:1440
components <- 1:3
time_axis <- scale_x_continuous(breaks = c(1, 481, 721, 1021, 1201), labels = c("12 AM", "8 AM", "12 PM", "5 PM", "8 PM"))
curve_colours <- c("Mean" = "grey35", "+ scorers" = "#B74F27", "- scorers" = "#2360A0")
curve_types <- c("Mean" = "solid", "+ scorers" = "dashed", "- scorers" = "dotted")

#### 2 Component effects: model-based and observed ####
# FPCA describes each participant as mean + score₁ × FPC1 + score₂ × FPC2 + …
# Row A: mean +/- 2 SD x eigenfunction (residual scale), the pattern each component adds or subtracts
# Row B: observed average day of participants with score >= +1 SD or <= -1 SD (asinh scale),  how the actual days of extreme scorers on that component compare with a typical day.
component_effects <- function(method, scores, efunctions, mu, variance) {
  labels <- paste0("FPC", components, " (", sprintf("%.1f", 100 * variance[components]), "%)")
  model_curves <- list()
  group_curves <- list()
  group_activity <- list()
  for (j in seq_along(components)) {
    k <- components[j]
    score_sd <- sd(scores[, k])
    model_curves[[j]] <- data.frame(component = labels[j],
      minute = time, mean = mu, 
      high = mu + 2 * score_sd * efunctions[, k],
      low = mu - 2 * score_sd * efunctions[, k])
    z <- as.numeric(scale(scores[, k]))
    high <- z >= 1
    low <- z <= -1
    group_curves[[j]] <- data.frame(component = labels[j], minute = time, mean = colMeans(subject_observed), 
    high = colMeans(subject_observed[high, , drop = FALSE]), low = colMeans(subject_observed[low, , drop = FALSE]))
    daily_high <- rowMeans(subject_observed[high, , drop = FALSE])
    daily_low <- rowMeans(subject_observed[low, , drop = FALSE])
    group_activity[[j]] <- data.frame(dataset, cosinor_fit, method, component = paste0("FPC", k),
    n_high = sum(high), n_low = sum(low),
    mean_activity_all = mean(rowMeans(subject_observed)),
    mean_activity_high = mean(daily_high),
    mean_activity_low = mean(daily_low),
    welch_p = t.test(daily_high, daily_low)$p.value)
  }
  long <- bind_rows(
    bind_rows(model_curves) %>% mutate(row = "A. Model (asinh residual)"),
    bind_rows(group_curves) %>% mutate(row = paste0("B. Observed (", unit, ")"))
  ) %>%
    pivot_longer(c(mean, high, low), names_to = "curve") %>%
    mutate(curve = factor(c(mean = "Mean", high = "+ scorers", low = "- scorers")[curve], names(curve_colours)), component = factor(component, labels),row = factor(row, unique(row)))
  p <- ggplot(long, aes(minute, value, colour = curve, linetype = curve)) +
    geom_hline(data = data.frame(row = levels(long$row)[1], yintercept = 0), aes(yintercept = yintercept), colour = "grey80") +
    geom_line(linewidth = 0.6) +
    facet_grid(row ~ component, scales = "free_y", switch = "y") +
    scale_colour_manual(values = curve_colours) +
    scale_linetype_manual(values = curve_types) +
    time_axis + theme_bw() +
    theme(legend.position = "bottom", strip.placement = "outside", strip.text.y.left = element_text(angle = 90)) +
    labs(title = paste0(toupper(dataset), ": ", method, " components of cosinor residuals"), subtitle = paste("A: mean +/- 2 SD x component.", "B: average observed day of participants scoring >= +1 SD (+) or <= -1 SD (-)."),
  x = "", y = NULL, colour = NULL, linetype = NULL)
  list(plot = p, activity = bind_rows(group_activity))
}

fpca_effects <- component_effects("FPCA", fpca_fit$scores, fpca_fit$efunctions, fpca_fit$mu, fpca_variance)
mfpca_effects <- component_effects("MFPCA level 1", mfpca_fit$scores$level1, mfpca_fit$efunctions$level1,  mfpca_fit$mu, mfpca_variance)
fpca_effects$plot
mfpca_effects$plot

#### 3 MFPCA level-2 (day-to-day) components ####
# Mean +/- 2 SD of daily scores
# A person's within-person SD on a component measures how widely their days spread between these curves.
level2_components <- 1:3
level2_share <- mfpca_fit$evalues$level2 / sum(mfpca_fit$evalues$level2)
plot_data <- bind_rows(lapply(level2_components, function(k) {
  score_sd <- sd(mfpca_fit$scores$level2[, k])
  data.frame(component = paste0("Level-2 FPC", k, " (", sprintf("%.1f", 100 * level2_share[k]), "% of retained level-2 variance)"),
  minute = time, mean = mfpca_fit$mu,
  plus = mfpca_fit$mu + 2 * score_sd * mfpca_fit$efunctions$level2[, k],
  minus = mfpca_fit$mu - 2 * score_sd * mfpca_fit$efunctions$level2[, k])
})) %>%
  mutate(component = factor(component, unique(component))) %>%
  pivot_longer(c(minus, mean, plus), names_to = "curve") %>%
  mutate(curve = factor(curve, c("minus", "mean", "plus"),
                        c("Mean - 2 SD daily score", "Mean", "Mean + 2 SD daily score")))
p <- ggplot(plot_data, aes(minute, value, colour = curve)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ component, ncol = 1) +
  scale_colour_manual(values = c("#2360A0", "grey30", "#B74F27")) +
  time_axis + theme_bw() +
  theme(legend.position = "bottom", strip.text = element_text(hjust = 0)) +
  labs(title = paste0(toupper(dataset), ": MFPCA level-2 components"),
       x = "", y = "Cosinor residual (asinh units)", colour = NULL)
p
# eigenfunction signs are arbitrary, so only magnitude and relative direction within a component are interpretable.
