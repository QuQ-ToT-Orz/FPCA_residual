#### Comparison with standard FPCA ####
#### 0 Packages ####
gc()
library(dplyr)
library(refund)
library(ggplot2)

#### 1 Settings ####
dir_path <- "./results/comparison/"
dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
load("./results/clustering/clustering.RData")       
load("./results/mortality/cohorts_weighted.RData")     
time <- 1:1440
nf <- function(x) x / sqrt(sum(x^2))
cs <- function(a, b) abs(sum(nf(a) * nf(b)))

#### 2 Harmonic ####
basis <- function(h) cbind(cos(2 * pi * h * (time - 1) / 1440), sin(2 * pi * h * (time - 1) / 1440))
# harmonic structure and component correspondence 
# Share of a (centred) eigenfunction at each period: 24, 12, 8, 6, 4.8, 4 hours
period_share <- function(phi, H = 1:6) {
  phi <- phi - mean(phi)
  e <- sapply(H, function(h) { X <- basis(h); sum((X %*% qr.solve(X, phi))^2) })
  setNames(100 * e / sum(phi^2), paste0(24 / H, "h"))
}
participant_mean <- function(m, id) rowsum(m, id) / as.vector(table(id))

harmonics <- list(); content <- list(); correspondence <- list(); 
cosinor_share <- list()
for (ds in c("count", "mims")) {
  suffix <- if (ds == "mims") "new" else "old"
  e <- new.env()
  load(paste0("./data/", ds, "/Act_Analysis_asinh_", suffix, ".RData"), envir = e)
  load(paste0("./results/cosinor_residual_fpca/", ds, "/residual_fpca_weekly_", suffix, ".RData"), envir = e)
  a <- e$Act_Analysis_asinh
  R <- e$subject_residual
  Y <- participant_mean(as.matrix(a[, paste0("MIN", time)]), as.character(a$SEQN))[rownames(R), ]

  #### 1 Harmonic coefficients of each participant's average day ####
  co <- do.call(cbind, lapply(1:4, function(h) Y %*% basis(h) * 2 / 1440))
  colnames(co) <- paste0(rep(c("c", "s"), 4), rep(1:4, each = 2))
  harmonics[[ds]] <- data.frame(
    dataset = ds, period = c("24h", "12h", "8h", "6h"),
    variance = sapply(1:4, function(h) var(co[, 2 * h - 1]) + var(co[, 2 * h])),
    # between-subject harmonic variability
    # 24-hour, 12-hour, 8-hour, or 6-hour harmonic contribution to curve variance
    variance_cos = sapply(1:4, function(h) var(co[, 2 * h - 1])),
    variance_sin = sapply(1:4, function(h) var(co[, 2 * h])),
    cancor_with_24h = c(NA, sapply(2:4, function(h) cancor(co[, c("c1", "s1")], co[, c(2 * h - 1, 2 * h)])$cor[1])))
    # variations in the 24-hour harmonic coefficients correlated with variations in the 12-hour, 8-hour, or 6-hour harmonic coefficients across participants

  #### 2 Harmonic content of standard FPC1-FPC6 and residual FPC1-FPC3 ####
  content[[ds]] <- bind_rows(
    lapply(1:6, function(k) data.frame(dataset = ds, method = "standard", component = k, t(period_share(fpca[[ds]]$efunctions[, k])))),
    lapply(1:3, function(k) data.frame(dataset = ds, method = "residual", component = k, t(period_share(cohorts[[ds]]$efunctions[, k])))))

  #### 3 Full-cohort correspondence: residual FPC1-FPC3 vs standard FPC1-FPC13 ####
  correspondence[[ds]] <- bind_rows(lapply(1:3, function(k) {
    s <- sapply(1:13, function(j) cs(cohorts[[ds]]$efunctions[, k], fpca[[ds]]$efunctions[, j])) # cosine similarity
    o <- order(-abs(s))[1:2]
    q <- qr.Q(qr(fpca[[ds]]$efunctions[, o]))
    data.frame(dataset = ds, residual = k, best_standard = o[1], similarity = s[o[1]],
               second_standard = o[2], similarity_second = s[o[2]],
               similarity_pair = sqrt(sum(crossprod(q, nf(cohorts[[ds]]$efunctions[, k]))^2)),
               standard_variance = 100 * fpca[[ds]]$variance[o[1]])
  }))

  #### 4 Share of the standard FPC4 score ####
  # Average day = cosinor part + residual; project both parts onto the standard FPC4 eigenfunction
  phi <- nf(fpca[[ds]]$efunctions[, 4])
  pc <- as.vector((Y - R) %*% phi); pr <- as.vector(R %*% phi)
  cosinor_share[[ds]] <- data.frame(dataset = ds,
    share_cosinor = var(pc) / var(pc + pr))
    # cor_score_residual_part = cor(fpca[[ds]]$scores[rownames(R), 4], pr)
  rm(a, e, R, Y, co); gc()
}
harmonics <- bind_rows(harmonics)
content <- bind_rows(content)
correspondence <- bind_rows(correspondence)
cosinor_share <- bind_rows(cosinor_share)
print(harmonics %>% mutate(across(where(is.numeric), ~ round(.x, 3))), row.names = FALSE)
print(content %>% mutate(across(where(is.numeric), ~ round(.x))), row.names = FALSE)
print(correspondence %>% mutate(across(where(is.double), ~ round(.x, 2))), row.names = FALSE)
print(cosinor_share %>% mutate(across(where(is.numeric), ~ round(.x, 2))), row.names = FALSE)

#### 3 Stability of the components ####
# FPC5, FPC4, FPC6 for residual FPC1, FPC2, FPC3
# Similarity: absolute inner product of the normalised eigenfunctions
sizes <- c(300, 500, 1000, 2000)
n_samples <- 30
standard_rank <- c(5, 4, 6)   # standard component matching residual FPC1, FPC2, FPC3
compare <- function(Y, R, target, idx) {
  fr <- fpca.face(Y = Y[idx, ], argvals = time)
  fs <- fpca.face(Y = R[idx, ], argvals = time)
  do.call(rbind, lapply(1:3, function(k) {
    rs <- sapply(1:5, function(j) cs(fs$efunctions[, j], target[, k]))
    ss <- sapply(1:10, function(j) cs(fr$efunctions[, j], target[, k]))
    data.frame(component = k, res_same = rs[k], res_best_rank = which.max(rs),
               std_fixed = ss[standard_rank[k]], std_best = max(ss), std_best_rank = which.max(ss))
  }))
}

## Subsamples and waves ##
subsamples <- list(); waves <- list()
for (ds in c("count", "mims")) {
  suffix <- if (ds == "mims") "new" else "old"
  e <- new.env()
  load(paste0("./data/", ds, "/Act_Analysis_asinh_", suffix, ".RData"), envir = e)
  load(paste0("./results/cosinor_residual_fpca/", ds, "/residual_fpca_weekly_", suffix, ".RData"), envir = e)
  a <- e$Act_Analysis_asinh
  R <- e$subject_residual
  Y <- participant_mean(as.matrix(a[, paste0("MIN", time)]), as.character(a$SEQN))[rownames(R), ]
  target <- e$fpca_fit$efunctions[, 1:3]
  wave <- cohorts[[ds]]$d$SDDSRVYR[match(rownames(R), cohorts[[ds]]$d$SEQN)]
  rm(a, e); gc()
  for (n in sizes) for (b in seq_len(n_samples)) {
    set.seed(1000*n + b)
    subsamples[[length(subsamples) + 1]] <- cbind(dataset = ds, n = n, sample = b, compare(Y, R, target, sample(nrow(R), n)))
  }
  for (w in sort(unique(wave))) {
    waves[[length(waves) + 1]] <- cbind(dataset = ds, wave = w, n = sum(wave == w), compare(Y, R, target, which(wave == w)))
  }
  rm(Y, R); gc()
}
subsamples <- do.call(rbind, subsamples)
waves <- do.call(rbind, waves)

# Each point is one random sample; bars are medians; dashed line: similarity 0.9.
cohort_label <- c(count = "Hip counts, 2003-2006", mims = "Wrist MIMS, 2011-2014")
pair_label <- c("Residual FPC1 vs standard FPC5", "Residual FPC2 vs standard FPC4", "Residual FPC3 vs standard FPC6")
colours <- c("#2a78d6", "#eb6834")

stability_plot <- function(d, y_label) {
  ggplot(d, aes(n, sim, colour = method)) +
    geom_hline(yintercept = 0.9, linetype = "dashed", colour = "#c3c2b7") +
    geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6, seed = 1), size = 1, alpha = 0.6) +
    stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.4, position = position_dodge(width = 0.6), show.legend = FALSE) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_colour_manual(values = colours, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "top", legend.justification = "left", panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold")) +
    labs(x = "Participants per sample (30 random samples each)", y = y_label)
}
long <- bind_rows(
  subsamples %>% transmute(dataset, n, component, method = "Residual FPCA", sim = res_same),
  subsamples %>% transmute(dataset, n, component, method = "Standard FPCA (component at its full-cohort rank)", sim = std_fixed)) %>%
  mutate(cohort = factor(cohort_label[dataset], cohort_label), n = factor(n), pair = factor(pair_label[component], pair_label),
         method = factor(method, c("Residual FPCA", "Standard FPCA (component at its full-cohort rank)")))

p <- stability_plot(filter(long, component == 2) %>%
                           mutate(method = factor(ifelse(method == "Residual FPCA", "Residual FPCA: FPC2", "Standard FPCA: FPC4"),
                                                  c("Residual FPCA: FPC2", "Standard FPCA: FPC4"))),
                         "Similarity to the full-cohort FPC2") +
  facet_wrap(~ cohort)
p

p <- stability_plot(long, "Similarity to the full-cohort component") + facet_grid(cohort ~ pair)
p

#### 4 Example: one n = 300 sample per cohort (lower quartile of standard FPC4 similarity) ####
hour <- (time - 1) / 60
# Sample at the lower quartile of standard FPC4 similarity
pick <- subsamples %>% filter(n == 300, component == 2) %>% group_by(dataset) %>%
  slice(which.min(abs(std_fixed - quantile(std_fixed, 0.25)))) %>% ungroup()

curves <- list()
for (i in seq_len(nrow(pick))) {
  ds <- pick$dataset[i]
  suffix <- if (ds == "mims") "new" else "old"
  e <- new.env()
  load(paste0("./data/", ds, "/Act_Analysis_asinh_", suffix, ".RData"), envir = e)
  load(paste0("./results/cosinor_residual_fpca/", ds, "/residual_fpca_weekly_", suffix, ".RData"), envir = e)
  a <- e$Act_Analysis_asinh
  R <- e$subject_residual
  Y <- participant_mean(as.matrix(a[, paste0("MIN", time)]), as.character(a$SEQN))[rownames(R), ]
  target <- nf(e$fpca_fit$efunctions[, 2])
  set.seed(1000 * 300 + pick$sample[i]); idx <- sample(nrow(R), 300)  
  fr <- fpca.face(Y = Y[idx, ], argvals = time)                    
  fs <- fpca.face(Y = R[idx, ], argvals = time)

  r2 <- nf(fs$efunctions[, 2]); r2 <- r2 * sign(sum(r2 * target))
  p4 <- nf(fr$efunctions[, standard_rank[2]]); p4 <- p4 * sign(sum(p4 * target))
  curves[[ds]] <- data.frame(
    cohort = cohort_label[[ds]], hour = rep(hour, 3), type = rep(c("ref", "res", "raw"), each = 1440),
    curve = rep(c("Full-cohort pattern",
                  sprintf("Residual FPC2, n = 300 (%.2f)", cs(r2, target)),
                  sprintf("Standard FPC4, n = 300 (%.2f)", cs(p4, target))), each = 1440),
    phi = c(target, r2, p4))
  rm(a, e, Y, R); gc()
}
b <- bind_rows(curves) %>% mutate(cohort = factor(cohort, cohort_label))

p <- ggplot(b, aes(hour, phi, group = curve, colour = type, linetype = type)) +
  geom_hline(yintercept = 0, colour = "#c3c2b7") +
  geom_line(linewidth = 0.6) +
  geom_text(data = b %>% filter(type != "ref") %>% distinct(cohort, curve, type) %>%
              mutate(hour = 0.3, phi = ifelse(type == "res", 0.075, 0.064)),
            aes(label = curve), hjust = 0, size = 2.9, show.legend = FALSE) +
  facet_wrap(~ cohort) +
  scale_colour_manual(values = c(ref = "#898781", res = colours[1], raw = colours[2]), guide = "none") +
  scale_linetype_manual(values = c(ref = "dashed", res = "solid", raw = "solid"), guide = "none") +
  scale_x_continuous(breaks = seq(0, 24, 4), labels = c("12 AM", "4 AM", "8 AM", "12 PM", "4 PM", "8 PM", "12 AM")) +
  scale_y_continuous(limits = c(-0.075, 0.08)) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", hjust = 0)) +
  labs(x = NULL, y = "Eigenfunction (unit norm)")
p
