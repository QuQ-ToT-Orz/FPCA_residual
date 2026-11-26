#### 0 Packages ####
gc()
library(dplyr)
library(ggplot2)
library(cluster)
library(refund)

#### 1 Settings ####
dir_path <- "./results/clustering/"
dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
load("./results/mortality/cohorts_weighted.RData")
k_max <- 8
n_reference <- 50   # uniform reference sets for the gap statistic
seed <- 123
variance_target <- 0.90
time <- 1:1440
nf <- function(x) x / sqrt(sum(x^2))
fit_kmeans <- function(x, k, nstart) kmeans(x, centers = k, nstart = nstart, iter.max = 1000, algorithm = "Lloyd")
converged <- function(km) (is.null(km$ifault) || km$ifault == 0) && km$iter < 1000
participant_mean <- function(m, id) rowsum(m, id) / as.vector(table(id))
cohort_label <- c(count = "Hip counts, 2003-2006", mims = "Wrist MIMS, 2011-2014")
cluster_colours <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442", "#999999")
hour_breaks <- c(0, 4, 8, 12, 16, 20, 24)
hour_labels <- c("12 AM", "4 AM", "8 AM", "12 PM", "4 PM", "8 PM", "12 AM")

#### 2 FPCA of participant-average asinh activity ####
fpca <- list()
curves <- list()
for (ds in names(cohorts)) {
suffix <- if (ds == "mims") "new" else "old"
id <- cohorts[[ds]]$d$SEQN
load(paste0("./data/", ds, "/Act_Analysis_asinh_", suffix, ".RData"))
asinh_curve <- participant_mean(as.matrix(Act_Analysis_asinh[, paste0("MIN", time)]), as.character(Act_Analysis_asinh$SEQN))[id, ]

load(paste0("./data/", ds, "/Act_Analysis_", suffix, ".RData"))
raw <- as.matrix(Act_Analysis[, paste0("MIN", time)])
raw[is.na(raw)] <- 0
raw_curve <- participant_mean(raw, as.character(Act_Analysis$SEQN))[id, ]

f <- fpca.face(Y = asinh_curve, argvals = time, pve = variance_target)
curves[[ds]] <- list(asinh = asinh_curve, raw = raw_curve, model = f$Yhat)
# Eigenfunction signs are arbitrary: make each one's largest absolute value positive
flip <- sign(apply(f$efunctions, 2, function(phi) phi[which.max(abs(phi))]))
f$efunctions <- sweep(f$efunctions, 2, flip, "*")
f$scores <- sweep(f$scores, 2, flip, "*")
rownames(f$scores) <- id
colnames(f$scores) <- paste0("FPC", seq_len(ncol(f$scores)))
# Each score is exactly shrinkage x projection; mu and shrinkage score new participants (Shiny_tool)
projection <- sweep(asinh_curve, 2, f$mu) %*% f$efunctions
fpca[[ds]] <- list(efunctions = f$efunctions, scores = f$scores, mu = f$mu, 
  shrinkage = colSums(f$scores * projection) / colSums(projection^2),
  variance = f$evalues / sum(fpca.face(Y = asinh_curve, argvals = time)$evalues))
rm(asinh_curve, raw_curve, f)
  gc()
}
# Align wrist FPC1-FPC3 signs to hip, so the FPC plots are comparable (k-means does not depend on signs)
s <- sign(sapply(1:3, function(k) sum(nf(fpca$count$efunctions[, k]) * nf(fpca$mims$efunctions[, k]))))
fpca$mims$efunctions[, 1:3] <- sweep(fpca$mims$efunctions[, 1:3], 2, s, "*")
fpca$mims$scores[, 1:3] <- sweep(fpca$mims$scores[, 1:3], 2, s, "*")

#### 3 K-means: K by the gap statistic (Tibs2001SEmax) ####
# Clustering uses unscaled functional scores, preserving relative component variances.  
clustering <- list()
k_diagnostics <- list()
cluster_summary <- list()
for (ds in names(fpca)) {
  x <- fpca[[ds]]$scores
  set.seed(seed)
  gap <- clusGap(x, function(x, k) fit_kmeans(x, k, 10), K.max = k_max, B = n_reference, d.power = 2, verbose = FALSE)
  best_k <- maxSE(gap$Tab[, "gap"], gap$Tab[, "SE.sim"], "Tibs2001SEmax")
  if (best_k == 1) stop(ds, ": gap statistic selects K = 1 (no cluster structure)")
  # K-means for every K, with convergence checks; silhouette for each
  fits <- lapply(2:k_max, function(k) { set.seed(seed); fit_kmeans(x, k, 50) })
  x_dist <- dist(x)
  k_diagnostics[[ds]] <- data.frame(
    dataset = ds, k = 1:k_max, gap = gap$Tab[, "gap"], gap_se = gap$Tab[, "SE.sim"],
    silhouette = c(NA, sapply(fits, function(f) mean(silhouette(f$cluster, x_dist)[, "sil_width"]))),
    # Convergence is a property of each k-means fit (no k-means fit for K = 1)
    converged = c(NA, ifelse(sapply(fits, converged), "yes", "no")),
    selected = 1:k_max == best_k)
  km <- fits[[best_k - 1]]
  # Label clusters by decreasing size so cluster 1 is the largest
  cluster <- match(km$cluster, order(-km$size))
  clustering[[ds]] <- list(k = best_k, gap = gap, cluster = setNames(cluster, rownames(x)), centers = km$centers[order(-km$size), , drop = FALSE])
  cluster_summary[[ds]] <- data.frame(
    dataset = ds, n_participants = nrow(x), n_components = ncol(x),
    variance_explained = sum(fpca[[ds]]$variance[1:ncol(x)]), 
    selected_k = best_k,
    silhouette = k_diagnostics[[ds]]$silhouette[best_k], converged = k_diagnostics[[ds]]$converged[best_k],
    cluster_sizes = paste(tabulate(cluster), collapse = "/"))
  rm(x_dist, fits)
}
k_diagnostics <- bind_rows(k_diagnostics)
cluster_summary <- bind_rows(cluster_summary)
print(k_diagnostics %>% mutate(across(where(is.double), ~ round(.x, 3))), row.names = FALSE)
print(cluster_summary %>% mutate(across(where(is.double), ~ round(.x, 3))), row.names = FALSE)
# dataset n_participants n_components variance_explained selected_k silhouette converged   cluster_sizes
#    count           2615           13              0.916          3      0.198       yes    1050/863/702
#     mims           3208           13              0.918          4      0.176       yes 958/872/847/531

#### 4 Clustering diagnostics: gap statistic and silhouette ####
panel <- function(ds) paste0(cohort_label[ds], " (K = ", clustering[[ds]]$k, ")")
diag_data <- k_diagnostics %>% mutate(cohort = factor(cohort_label[dataset], cohort_label))
p_gap <- ggplot(diag_data, aes(k, gap)) +
  geom_line(colour = "grey45") +
  geom_errorbar(aes(ymin = gap - gap_se, ymax = gap + gap_se), width = 0.25, colour = "grey45") +
  geom_point(colour = "grey25", size = 1.8) +
  geom_point(data = filter(diag_data, selected), colour = "#B74F27", size = 3.5) +
  geom_text(data = filter(diag_data, selected), aes(label = paste0("K = ", k)), colour = "#B74F27", vjust = -1.2, size = 3.5) +
  facet_wrap(~ cohort, scales = "free_y") +
  scale_x_continuous(breaks = 1:k_max) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.caption = element_text(hjust = 0.5)) +
  labs(x = "Number of clusters K", y = "Gap statistic",
       caption = "Bars: gap +/- 1 SE (50 uniform reference sets). \nHighlighted: Tibshirani (2001) 1-SE rule.")
p_gap

# Average silhouette width for K = 2-8 (undefined for K = 1)
p_sil <- ggplot(filter(diag_data, k > 1), aes(k, silhouette)) +
  geom_line(colour = "grey45") +
  geom_point(colour = "grey25", size = 1.8) +
  geom_point(data = filter(diag_data, selected), colour = "#B74F27", size = 3.5) +
  geom_text(data = filter(diag_data, selected), aes(label = paste0("K = ", k)), colour = "#B74F27", vjust = -1.2, size = 3.5) +
  facet_wrap(~ cohort) +
  scale_x_continuous(breaks = 2:k_max) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15))) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.caption = element_text(hjust = 0.5)) +
  labs(x = "Number of clusters K", y = "Average silhouette width",
       caption = "Highlighted: K chosen by the gap statistic. \nSilhouette below 0.25 indicates no substantial cluster structure.")
p_sil

#### 5 Clusters on FPC1 vs FPC2 ####
score_data <- bind_rows(lapply(names(clustering), function(ds)
  data.frame(panel = panel(ds), FPC1 = fpca[[ds]]$scores[, "FPC1"], FPC2 = fpca[[ds]]$scores[, "FPC2"],
             cluster = factor(clustering[[ds]]$cluster))))
centres <- score_data %>% group_by(panel, cluster) %>% summarise(FPC1 = mean(FPC1), FPC2 = mean(FPC2), .groups = "drop")
p_scores <- ggplot(score_data, aes(FPC1, FPC2, colour = cluster)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(size = 0.9, alpha = 0.45) +
  geom_point(data = centres, shape = 21, fill = "white", size = 4, stroke = 1.3) +
  geom_text(data = centres, aes(label = cluster), colour = "black", size = 2.8) +
  facet_wrap(~ panel, scales = "free") +
  scale_colour_manual(values = cluster_colours) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.caption = element_text(hjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(x = "FPC1 score", y = "FPC2 score", colour = "Cluster")
p_scores

#### 6 Cluster mean curves: model-fitted and observed ####
# Model: mean of the cluster's FPCA-fitted curves (asinh scale only).
# Observed: mean of the cluster's participant-average curves, on the asinh and raw scales.
profile_curves <- bind_rows(lapply(names(clustering), function(ds) {
  cl <- clustering[[ds]]$cluster
  cluster_means <- function(m, scale, curve) {
    means <- sapply(seq_len(clustering[[ds]]$k), function(j) colMeans(m[cl == j, , drop = FALSE]))
    data.frame(dataset = ds, scale = scale, curve = curve, cluster = rep(seq_len(ncol(means)), each = length(time)),
               hour = (time - 1) / 60, value = c(means))
  }
  bind_rows(cluster_means(curves[[ds]]$model, "asinh", "Model"), cluster_means(curves[[ds]]$asinh, "asinh", "Observed"), cluster_means(curves[[ds]]$raw, "raw", "Observed"))
}))
p_curves <- profile_curves %>%
  mutate(panel = paste0(cohort_label[dataset], "\n", curve, " (", scale, ")"), cluster = factor(cluster)) %>%
  ggplot(aes(hour, value, colour = cluster)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ panel, ncol = 3, scales = "free_y") +
  scale_x_continuous(breaks = hour_breaks, labels = hour_labels) +
  scale_colour_manual(values = cluster_colours) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
        plot.caption = element_text(hjust = 0.5)) +
  labs(x = NULL, y = "Mean activity per minute", colour = "Cluster")
p_curves

#### 7 Save ####
assignments <- bind_rows(lapply(names(clustering), function(ds)
  data.frame(dataset = ds, SEQN = names(clustering[[ds]]$cluster), cluster = clustering[[ds]]$cluster)))
save(fpca, clustering, k_diagnostics, cluster_summary, assignments, profile_curves, file = paste0(dir_path, "clustering.RData"))
