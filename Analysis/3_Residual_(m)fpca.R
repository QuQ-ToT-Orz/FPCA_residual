#### 0 Packages ####
gc()
library(dplyr)
library(tidyr)
library(ggplot2)
library(refund)

#### 1 Load the data ####
dataset <- "mims"
dataset <- "count"
cosinor_fit <- "weekly"   # weekly (default) or "daily"
expanded <- TRUE   # FALSE for the default run
out_suffix <- if (expanded) "_expanded" else ""
suffix <- paste0(if (dataset == "mims") "new" else "old", out_suffix)
load(paste0("./data/", dataset, "/Act_Analysis_asinh_", suffix, ".RData"))
Act_Analysis_asinh <- Act_Analysis_asinh %>% arrange(SEQN, WEEKDAY)
observed <- as.matrix(Act_Analysis_asinh[, paste0("MIN", 1:1440)])
index <- Act_Analysis_asinh[, c("SEQN", "WEEKDAY")]

set.seed(123)
sample_id <- Act_Analysis_asinh %>%
  group_by(Gender, Race, MobilityProblem) %>%
  summarise(ID = list(unique(SEQN)), .groups = "keep") %>%
  mutate(sampled = lapply(ID, function(x) sample(x, size = round(length(x) * 0.15)))) %>%
  unnest(cols = c(sampled))
id <- c(head(sample_id, n = 1)$sampled, tail(sample_id, n = 1)$sampled)
id
rm(sample_id)

time <- 1:1440
day_id <- as.character(index$SEQN)
subject_id <- sort(unique(day_id))
variance_target <- 0.90

dir_path <- paste0("./results/cosinor_residual_fpca/", dataset, "/")
dir.create(dir_path, recursive = TRUE, showWarnings = TRUE)

#### 2 Cosinor residuals ####
load(paste0("./results/cosinor_analysis/", dataset, "/daily_cosinor_parameters_", suffix, ".RData"))
load(paste0("./results/cosinor_analysis/", dataset, "/weekly_cosinor_parameters_", suffix, ".RData"))
cosinor_curves <- function(p) {
p$mes + p$amp * cos(2 * pi * outer(p$acrotime, time, function(a, t) t - a) / 1440)
}
if (cosinor_fit == "daily") {
  p <- daily_parameters$cosinor
} else {
  p <- weekly_parameters$cosinor[match(day_id, as.character(weekly_parameters$cosinor$SEQN)), ]
}
# Per day: activity minus that day's curve (daily) or the participant's weekly curve (weekly)
cosinor_residual <- observed - cosinor_curves(p)
rm(observed, p)
# Per participant: average of the 7 daily residuals. The same for daily and weekly fits,
# Daily fit: average over the 7 days of (activity − that day's curve) = (average activity) − (average of the 7 daily curves)
# Weekly fit: average over the 7 days of (activity − the weekly curve) = (average activity) − (weekly curve)
subject_residual <- t(sapply(split(seq_len(nrow(cosinor_residual)), day_id),
                             function(inx) colMeans(cosinor_residual[inx, , drop = FALSE])))

#### 3 FPCA of participant-average residuals ####
# Internal FACE covariance smoothing with 35 knots and data-selected penalties
# Default B-spline boundaries, no periodic constraint.
fpca_fit <- fpca.face(Y = subject_residual, argvals = time, pve = variance_target)
# Eigenfunction signs are arbitrary. 
# Participant-average residuals sum to 0 (each day does under the daily fit;
# the 7 days together under the weekly fit), so the eigenfunctions integrate to 0 too; 
# instead make each one's largest absolute value positive, and flip scores to match.
flip <- sign(apply(fpca_fit$efunctions, 2, function(f) f[which.max(abs(f))]))
fpca_fit$efunctions <- sweep(fpca_fit$efunctions, 2, flip, "*")
fpca_fit$scores <- sweep(fpca_fit$scores, 2, flip, "*")
rownames(fpca_fit$scores) <- subject_id
# fpca.face does not return its shrinkage, evalues / (evalues + sigma2). 
# Each score is exactly shrinkage x projection
projection <- sweep(subject_residual, 2, fpca_fit$mu) %*% fpca_fit$efunctions
fpca_fit$shrinkage <- colSums(fpca_fit$scores * projection) / colSums(projection^2)
fpca_variance <-fpca_fit$pve * fpca_fit$evalues / sum(fpca_fit$evalues)  # approx. share of total variance

#### 4 Multilevel FPCA of daily residuals ####
# Level 1 = between participants (used for clustering); level 2 = day-to-day within participants.
mfpca_fit <- mfpca.face(
  Y = cosinor_residual,
  id = factor(day_id, levels = subject_id),
  visit = index$WEEKDAY,
  argvals = time, weight = "subj",
  pve = variance_target
)
# Match each level-1 eigenfunction's sign to the corresponding FPCA eigenfunction.
k <- min(ncol(mfpca_fit$efunctions$level1), ncol(fpca_fit$efunctions))
flip <- c(sign(colSums(mfpca_fit$efunctions$level1[, 1:k] * fpca_fit$efunctions[, 1:k])),
          rep(1, ncol(mfpca_fit$efunctions$level1) - k))
mfpca_fit$efunctions$level1 <- sweep(mfpca_fit$efunctions$level1, 2, flip, "*")
mfpca_fit$scores$level1 <- sweep(mfpca_fit$scores$level1, 2, flip, "*")
rownames(mfpca_fit$scores$level1) <- subject_id
mfpca_variance <- mfpca_fit$pve * mfpca_fit$evalues$level1 / sum(mfpca_fit$evalues$level1)

#### 5 Fitted participant curves ####
# subject_residual is noisy minute by minute; these reconstructions show the pattern the FPCs represent. 
# MFPCA: average of each participant's level-1 day curves (mean + level-1 part).
fitted_subject_residual <- list(
  fpca = fpca_fit$Yhat,
  mfpca = sweep(mfpca_fit$scores$level1 %*%t(mfpca_fit$efunctions$level1), 2, mfpca_fit$mu, "+")
)
rownames(fitted_subject_residual$fpca) <- subject_id
rownames(fitted_subject_residual$mfpca) <- subject_id

#### 6 Plots ####
# eigenfunctions describe which residual temporal activity patterns vary across participants, 
# FPC scores describe how strongly each participant exhibits that pattern, 
# and eigenvalues measure how much variation that pattern explains across participants.
# First three eigenfunctions
# Panel label: method, FPC and its share of variance
fpc_label <- function(method, variance) paste0(method, ": FPC", 1:3, " (", sprintf("%.1f", 100 * variance[1:3]), "%)")
labels <- c(fpc_label("FPCA", fpca_variance), fpc_label("MFPCA level 1", mfpca_variance))
efunctions <- data.frame(
  label = factor(rep(labels, each = 1440), levels = labels),
  minute = time,
  value = c(fpca_fit$efunctions[, 1:3], mfpca_fit$efunctions$level1[, 1:3])
)

p <- ggplot(efunctions, aes(minute, value)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_line() + facet_wrap(~ label, nrow = 2, scales = "free_y") + theme_bw() +
  scale_x_continuous(breaks = c(1, 481, 721, 1021, 1201),
                     labels = c("12 AM", "8 AM", "12 PM", "5 PM", "8 PM")) +
  labs(x = "", y = "Eigenfunction")
p
# eigenfunction signs are arbitrary, so only magnitude and relative direction within a component are interpretable.

#### 7 Save ####
save(subject_residual, fitted_subject_residual, fpca_fit, fpca_variance, mfpca_fit, mfpca_variance, file = paste0(dir_path, "residual_fpca_", cosinor_fit, "_", suffix, ".RData"))

example_ids <- c(count = "30858", mims = "69831")
cohort_label <- c(count = "Hip counts", mims = "Wrist MIMS")
example_data <- function(ds) {
  sfx <- paste0(if (ds == "mims") "new" else "old", out_suffix)
  seqn <- example_ids[[ds]]
  participant <- paste0(cohort_label[[ds]], ", SEQN ", seqn)
  e <- new.env()
  load(paste0("./results/cosinor_residual_fpca/", ds, "/residual_fpca_", cosinor_fit, "_", sfx, ".RData"), envir = e)
  load(paste0("./data/", ds, "/Act_Analysis_asinh_", sfx, ".RData"), envir = e)
  load(paste0("./results/cosinor_analysis/", ds, "/", cosinor_fit, "_cosinor_parameters_", sfx, ".RData"), envir = e)
  a <- e$Act_Analysis_asinh %>% arrange(SEQN, WEEKDAY) 
  days <- which(as.character(a$SEQN) == seqn)
  curve_rows <- if (cosinor_fit == "daily") {
    e$daily_parameters$cosinor[days, ]
  } else {
    e$weekly_parameters$cosinor[rep(match(seqn, as.character(e$weekly_parameters$cosinor$SEQN)), length(days)), ]
  }
  list(
    average = data.frame(participant = participant, minute = time, residual = e$subject_residual[seqn, ], FPCA = e$fitted_subject_residual$fpca[seqn, ], MFPCA = e$fitted_subject_residual$mfpca[seqn, ]),
    daily = data.frame(participant = participant,
                       WEEKDAY = rep(as.character(a$WEEKDAY[days]), each = 1440),
                       minute = time,
                       residual = c(t(as.matrix(a[days, paste0("MIN", 1:1440)]) - cosinor_curves(curve_rows))))
  )
}
examples <- lapply(names(example_ids), example_data)
plot_data <- bind_rows(lapply(examples, `[[`, "average"))
daily_residuals <- bind_rows(lapply(examples, `[[`, "daily"))

# Average residual and its FPCA/MFPCA reconstructions
p <- ggplot(plot_data, aes(minute, residual)) +
  geom_line(colour = "grey65", linewidth = 0.2) +
  geom_line(aes(y = FPCA, colour = "FPCA")) +
  geom_line(aes(y = MFPCA, colour = "MFPCA level 1"), linetype = "dashed") +
  facet_wrap(~ participant, ncol = 1, scales = "free_y") + theme_bw() +   # hip and wrist are on different scales
  theme(legend.position = "bottom") +
  labs(x = "Minute of day", y = "Mean cosinor residual (asinh units)", colour = "")
p

# Daily residuals and their 7-day average, the curve analysed by FPCA
p <- ggplot(daily_residuals, aes(minute, residual)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_line(aes(colour = WEEKDAY), linewidth = 0.2, alpha = 0.7) +
  geom_line(data = plot_data, colour = "black", linewidth = 0.8) +
  facet_wrap(~ participant, ncol = 1, scales = "free_y") +
  scale_colour_brewer(palette = "Dark2", guide = "none") +
  scale_x_continuous(breaks = c(1, 481, 721, 1021, 1201),
                     labels = c("12 AM", "8 AM", "12 PM", "5 PM", "8 PM")) +
  theme_bw() +
  labs(x = "", y = "Cosinor residual (asinh units)",
       subtitle = "Thin coloured lines: individual days; thick black line: 7-day average")
p
