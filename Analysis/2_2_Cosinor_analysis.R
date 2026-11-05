#### 0 Packages ####
gc()
library(dplyr)
library(tidyr)
library(ggplot2)
library(ActCR)
# devtools::install_github("junruidi/ActCR")

#### 1 Load the data ####
dataset <- "mims" 
dataset <- "count"
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
id # 69831 30858

#### 2 Model functions ####
time <- 1:1440
weekdays <- c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")
fit_cosinor <- function(x) {
  fit <- ActCosinor(x, window = 1, export_ts = FALSE)
  p <- as.data.frame(fit$params)
  p$acrotime <- p$acrotime * 60
  p
}
cosinor_curve <- function(p) p$mes + p$amp * cos(2 * pi * (time - p$acrotime) / 1440)

# ActCR's exponential can overflow for very steep curves
fit_extended <- function(x) {
  fit <- ActExtendCosinor(x, window = 1, upper = c(Inf, Inf, 1, 300, 27), export_ts = FALSE)
  p <- as.data.frame(fit$params)
  p[c("acrotime", "UpMesor", "DownMesor")] <- p[c("acrotime", "UpMesor", "DownMesor")] * 60
  p
}
extended_curve <- function(p) {
  p$minimum + p$amp * plogis(p$beta * (cos(2 * pi * (time - p$acrotime) / 1440) - p$alpha))
}

dir_path <- paste0("./results/cosinor_analysis/", dataset, "/")
dir.create(dir_path, recursive = TRUE, showWarnings = TRUE)

#### 3 Daily fits: one cosinor per participant-day ####
empty <- matrix(NA_real_, nrow(observed), 1440, dimnames = dimnames(observed))
daily_fitted <- list(cosinor = empty, extended = empty)
daily_parameters <- list(cosinor = vector("list", nrow(observed)), extended = vector("list", nrow(observed)))
for (i in seq_len(nrow(observed))) {
  p <- fit_cosinor(observed[i, ])
  daily_fitted$cosinor[i, ] <- cosinor_curve(p)
  daily_parameters$cosinor[[i]] <- cbind(index[i, ], model = "cosinor", p)
  p <- fit_extended(observed[i, ])
  daily_fitted$extended[i, ] <- extended_curve(p)
  daily_parameters$extended[[i]] <- cbind(index[i, ], model = "extended", p)
  if (i %% 1000 == 0) message(dataset, " daily: fitted ", i, "/", nrow(observed))
}
daily_parameters <- lapply(daily_parameters, bind_rows)
save(daily_parameters, file = paste0(dir_path, "daily_cosinor_parameters_", suffix, ".RData"))

#### 4 Weekly fits: one cosinor per participant over all seven days ####
day_id <- as.character(index$SEQN)
subject_id <- sort(unique(day_id))
weekly_fitted <- list(cosinor = empty, extended = empty)
weekly_parameters <- list(cosinor = vector("list", length(subject_id)),
                          extended = vector("list", length(subject_id)))
for (i in seq_along(subject_id)) {
  inx <- which(day_id == subject_id[i])          # this participant's seven days
  x <- c(t(observed[inx, , drop = FALSE]))       # days concatenated into one series
  p <- fit_cosinor(x)
  weekly_fitted$cosinor[inx, ] <- rep(cosinor_curve(p), each = length(inx)) 
  weekly_parameters$cosinor[[i]] <- data.frame(SEQN = subject_id[i], model = "cosinor", p)
  p <- fit_extended(x)
  weekly_fitted$extended[inx, ] <- rep(extended_curve(p), each = length(inx))
  weekly_parameters$extended[[i]] <- data.frame(SEQN = subject_id[i], model = "extended", p)
  if (i %% 500 == 0) message(dataset, " weekly: fitted ", i, "/", length(subject_id))
}
weekly_parameters <- lapply(weekly_parameters, bind_rows)
save(weekly_parameters, file = paste0(dir_path, "weekly_cosinor_parameters_", suffix, ".RData"))

# ActCosinor is OLS with a fixed 24-hour period, 
# so the weekly fit is the cosinor of the participant's average day: 
# its curve, MESOR and cos/sin coefficients equal the average of the seven daily fits (amplitude and acrophase do not average). 
# Participant-average residuals are therefore the same under daily and weekly fits; 
# daily residuals differ, because day-to-day changes in MESOR, amplitude and timing stay in them.

#### 5 Plots ####
id # 69831 (new) 30858 (old)
curves <- function(parameters, rows, curve) t(sapply(rows, function(j) curve(parameters[j, ])))

# One example participant from each cohort
example_ids <- c(count = "30858", mims = "69831")
cohort_label <- c(count = "Hip counts", mims = "Wrist MIMS")
example_data <- function(ds) {
  sfx <- paste0(if (ds == "mims") "new" else "old", out_suffix)
  res_path <- paste0("./results/cosinor_analysis/", ds, "/")
  e <- new.env()
  load(paste0("./data/", ds, "/Act_Analysis_asinh_", sfx, ".RData"), envir = e)
  load(paste0(res_path, "daily_cosinor_parameters_", sfx, ".RData"), envir = e)
  load(paste0(res_path, "weekly_cosinor_parameters_", sfx, ".RData"), envir = e)
  a <- e$Act_Analysis_asinh %>% arrange(SEQN, WEEKDAY)   # same row order as the daily fits
  rows <- which(as.character(a$SEQN) == example_ids[[ds]])
  week <- match(example_ids[[ds]], e$weekly_parameters$cosinor$SEQN)
  data.frame(
    participant = paste0(cohort_label[[ds]], ", SEQN ", example_ids[[ds]]),
    WEEKDAY = factor(rep(as.character(a$WEEKDAY[rows]), each = 1440), weekdays),
    minute = time,
    observed = c(t(as.matrix(a[rows, paste0("MIN", 1:1440)]))),
    daily_cosinor = c(t(curves(e$daily_parameters$cosinor, rows, cosinor_curve))),
    daily_extended = c(t(curves(e$daily_parameters$extended, rows, extended_curve))),
    weekly_cosinor = cosinor_curve(e$weekly_parameters$cosinor[week, ]),
    weekly_extended = extended_curve(e$weekly_parameters$extended[week, ])
  )
}
plot_data <- bind_rows(lapply(names(example_ids), example_data))
p <- ggplot(plot_data, aes(minute, observed)) +
  geom_line(colour = "grey65", linewidth = 0.2) +
  # geom_line(aes(y = daily_cosinor, colour = "Daily cosinor")) +
  # geom_line(aes(y = daily_extended, colour = "Daily extended")) +
  geom_line(aes(y = weekly_cosinor, colour = "cosinor"), linetype = "dashed") +
  geom_line(aes(y = weekly_extended, colour = "cosinor extended"), linetype = "dashed") +
  facet_grid(participant ~ WEEKDAY, scales = "free_y") + theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Minute of day", y = "asinh(activity)", colour = "Model")
p

pdf(paste0("./results/cosinor_analysis/residual_diagnostics", out_suffix, ".pdf"), width = 8, height = 14)
par(mfrow = c(7, 2), oma = c(0, 0, 2, 0), mar = c(3.5, 3.5, 2, 1), mgp = c(2, 0.6, 0))
for (person in unique(plot_data$participant)) {
  for (day in levels(plot_data$WEEKDAY)) {
    d_day <- plot_data[plot_data$participant == person & plot_data$WEEKDAY == day, ]
    yhat <- d_day$weekly_cosinor
    r <- d_day$observed - yhat
    # plot(time, r, type = "l", xlab = "Minute of day", ylab = "Residual")
    # abline(h = 0, col = "red")
    # plot(yhat, r, pch = ".", xlab = "Fitted value", ylab = "Residual")
    # abline(h = 0, col = "red")
    qqnorm(r, pch = ".", main = paste(day, "- Normal Q-Q"))
    qqline(r, col = "red")
    if (sd(r) > 0) acf(r, lag.max = 180, main = "") else plot.new()
    title(paste(day, "- Residual ACF (lag in minutes)"), line = 0.8)
  }
  mtext(paste(person, "- cosinor"), outer = TRUE)
}
dev.off()
