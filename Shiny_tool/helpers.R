library(ActCR)
library(survey)

read_table <- function(path) {
  if (tolower(tools::file_ext(path)) == "rds") readRDS(path) else read.csv(path, check.names = FALSE)
}

# One participant's days x minutes matrix.
participant_matrix <- function(activity, id = NULL) {
  rows <- if ("SEQN" %in% names(activity)) as.character(activity$SEQN) == id else seq_len(nrow(activity))
  x <- as.matrix(activity[rows, paste0("MIN", 1:1440)])
  x[is.na(x) | x < 0] <- 0
  x
}

# Scores a participant against the fixed reference; no model is refitted.
score_participant <- function(x, ref, person) {
  y <- asinh(x)
  fit <- as.data.frame(ActCosinor(c(t(y)), window = 1, export_ts = FALSE)$params)
  acrophase <- fit$acrotime * 60
  cosinor <- fit$mes + fit$amp * cos(2 * pi * (1:1440 - acrophase) / 1440)
  residual <- colMeans(sweep(y, 2, cosinor))

  # FPC2 of the cosinor residuals
  fpc2 <- sum(ref$fpc2$efunction * (residual - ref$fpc2$mu)) * ref$fpc2$shrinkage
  fpc2_sd <- (fpc2 - ref$fpc2$mean) / ref$fpc2$sd

  # Nearest cluster, from the FPCA of average asinh activity
  curve <- colMeans(y)
  scores <- drop(crossprod(ref$cluster$efunctions, curve - ref$cluster$mu)) * ref$cluster$shrinkage
  distance <- sqrt(rowSums(sweep(ref$cluster$centers, 2, scores)^2))

  # Cardiovascular mortality; CI on the logit scale
  new <- data.frame(
    person, FPC2 = fpc2_sd,
    MESOR_z = (fit$mes - ref$cosinor[["mes_mean"]]) / ref$cosinor[["mes_sd"]],
    Amplitude_z = (fit$amp - ref$cosinor[["amp_mean"]]) / ref$cosinor[["amp_sd"]],
    cosinor_phase_sin = sin(2 * pi * acrophase / 1440), cosinor_phase_cos = cos(2 * pi * acrophase / 1440))
  p <- predict(ref$mortality$fit, newdata = new, type = "link", se.fit = TRUE)
  ci <- plogis(as.numeric(p) + c(0, -1.96, 1.96) * as.numeric(SE(p)))

  list(
    days = nrow(x),
    mesor = fit$mes, amplitude = fit$amp, acrophase = acrophase,
    fpc2 = fpc2_sd, fpc2_percentile = 100 * mean(ref$fpc2$scores <= fpc2),
    residual = residual, fpc2_curve = ref$fpc2$mu + fpc2 * ref$fpc2$efunction,
    cluster = which.min(distance), distance = distance, curve = curve,
    risk = ci, risk_percentile = 100 * mean(ref$mortality$risks <= ci[1]))
}

clock <- function(minute) sprintf("%02d:%02d", round(minute) %/% 60 %% 24, round(minute) %% 60)
hour_axis <- ggplot2::scale_x_continuous(breaks = seq(0, 24, 4), labels = c("12 AM", "4 AM", "8 AM", "12 PM", "4 PM", "8 PM", "12 AM"))
