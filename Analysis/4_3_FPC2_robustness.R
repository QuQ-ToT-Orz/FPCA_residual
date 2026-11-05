#### Robustness of the FPC2-mortality association ####
#### 0 Packages ####
gc()
library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(survey)

#### 1 Settings ####
dir_path <- "./results/mortality/"
load(paste0(dir_path, "cohorts_weighted.RData"))
base_terms <- c("Age", "Gender", "Race", "Income", "MESOR_z", "Amplitude_z", "cosinor_phase_sin", "cosinor_phase_cos")
activity_terms <- c("wear_z", "asinh_TAC", "MVPA_z")
fpc_terms <- c("FPC1", "FPC2", "FPC3")
main_terms <- c(base_terms, fpc_terms)
outcomes <- c(mortstat = "All-cause", cvd_death = "Cardiovascular", cancer_death = "Cancer")
methods <- c("unweighted", "weighted")
cohort_label <- c(count = "Hip counts, 2003-2006", mims = "Wrist MIMS, 2011-2014", pooled = "Pooled (fixed effect)")

ext_terms <- c("ext_minimum", "ext_alpha", "ext_log_beta")   # Extended cosinor beyond the ordinary cosinor
rhythm_terms <- c("IS", "IV")
lifestyle_terms <- c("BMI", "SmokeCigs", "DrinkStatus", "EducationAdult", "MobilityProblem")
disease_terms <- c("Diabetes", "CHF", "CHD", "Stroke",  "Cancer")

#### 2 correlation ####
# Spearman, with factors ranked by level order. Race is nominal, so it is left out; missing income is NA.
corr_vars <- c(setdiff(main_terms, "Race"), ext_terms, rhythm_terms, lifestyle_terms, activity_terms)
corr <- lapply(cohorts, function(x) {
  d <- x$d[corr_vars]
  d$Income[d$Income == "Missing"] <- NA
  d[] <- lapply(d, as.numeric)
  cor(d, method = "spearman", use = "pairwise.complete.obs")
})
high_corr <- bind_rows(lapply(names(corr), function(ds) {
  r <- corr[[ds]]
  i <- which(upper.tri(r) & abs(r) >= 0.7, arr.ind = TRUE)
  data.frame(dataset = ds, var1 = rownames(r)[i[, 1]], var2 = colnames(r)[i[, 2]], rho = r[i])
})) %>% arrange(dataset, desc(abs(rho)))
print(high_corr %>% mutate(rho = round(rho, 2)), row.names = FALSE)

corr_long <- bind_rows(lapply(names(corr), function(ds) {
  r <- corr[[ds]]
  i <- which(lower.tri(r), arr.ind = TRUE)
  data.frame(dataset = ds, var1 = rownames(r)[i[, 1]], var2 = colnames(r)[i[, 2]], rho = r[i])
})) %>%
  mutate(var1 = factor(var1, rev(corr_vars)), var2 = factor(var2, corr_vars),
         dataset = factor(cohort_label[dataset], cohort_label[names(corr)]))
p_corr <- ggplot(corr_long, aes(var2, var1, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(abs(rho) >= 0.7, sprintf("%.2f", rho), "")), size = 2, colour = "grey15") +
  facet_wrap(~dataset) +
  scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", midpoint = 0, limits = c(-1, 1),
                       name = "Spearman rho") +
  coord_fixed() +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
  labs(x = NULL, y = NULL)
p_corr

# Total activity and MVPA are highly correlated with each other (rho 0.88-0.97), and they and wear time overlap with the MESOR (rho up to 0.93), so each enters only as a separate adjustment, never together or in the combined model.
combined_terms <- c(ext_terms, rhythm_terms, "MobilityProblem")

specs <- list(
  list(group = "Adjustment", label = "Model 1", terms = main_terms),
  list(group = "Adjustment", label = "+ extended cosinor", terms = c(main_terms, ext_terms)),
  list(group = "Adjustment", label = "+ extended cosinor, IS and IV", terms = c(main_terms, ext_terms, rhythm_terms)),
  list(group = "Adjustment", label = "+ mobility", terms = c(main_terms, "MobilityProblem")),
  list(group = "Adjustment", label = "+ all above (Model 2)", terms = c(main_terms, combined_terms)),
  list(group = "Adjustment", label = "Model 2 + wear time", terms = c(main_terms, combined_terms, "wear_z")),
  # list(group = "Adjustment", label = "+ total activity", terms = c(main_terms, "asinh_TAC")),
  list(group = "Adjustment", label = "Model 2 + MVPA", terms = c(main_terms, combined_terms, "MVPA_z")),
  list(group = "Adjustment", label = "Model 2 + BMI, smoking, alcohol, education)",
  terms = c(main_terms, combined_terms, lifestyle_terms)),
  list(group = "Adjustment", label = "Model 2 + all above", terms = c(main_terms, combined_terms, lifestyle_terms, "wear_z", "MVPA_z")),
  list(group = "Sensitivity", label = "Deaths in first 2 years excluded", terms = main_terms, early = TRUE))
all_terms <- unique(unlist(lapply(specs, `[[`, "terms")))

#### 2 Model functions ####
designs <- lapply(cohorts, function(x) svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~w_ipw,
                                                 nest = TRUE, data = x$d))
cox <- function(d, terms, event, method, ds) {
  f <- as.formula(paste0("Surv(permth_exm / 12, ", event, ") ~ ", paste(terms, collapse = " + ")))
  if (method == "unweighted") return(coxph(f, data = d, ties = "efron"))
  svycoxph(f, design = subset(designs[[ds]], SEQN %in% d$SEQN))
}
estimate <- function(f, term) {
  b <- unname(coef(f)[term])
  se <- sqrt(unname(vcov(f)[term, term]))
  c(b = b, se = se, p = 2 * pnorm(-abs(b / se)))
}

#### 3 Fit every specification ####
robust <- list()
for (ev in names(outcomes)) for (ds in names(cohorts)) {
  d <- cohorts[[ds]]$d
  d <- d[complete.cases(d[, c("permth_exm", ev, all_terms)]), ]
  for (m in methods) {
    for (s in specs) {
      dd <- d
      if (isTRUE(s$early)) dd <- dd[!(dd$mortstat == 1 & dd$permth_exm < 24), ]
      tm <- if (is.null(s[["term"]])) "FPC2" else s[["term"]]   # [[ ]]: s$term would partially match s$terms
      f <- cox(dd, s$terms, ev, m, ds)
      k <- estimate(f, tm)
      robust[[length(robust) + 1]] <- data.frame(
        outcome = ev, dataset = ds, method = m, group = s$group, label = s$label, n = nrow(dd), deaths = sum(dd[[ev]]),
        b = k[["b"]], se = k[["se"]], p = k[["p"]],
        ph_p = if (m == "unweighted") cox.zph(f)$table[tm, "p"] else NA_real_) # proportional-hazards diagnostics
    }
  }
}
robust <- bind_rows(robust)

# Pooled across cohorts: fixed-effect inverse-variance; heterogeneity: Cochran's Q
pooled <- robust %>%
  group_by(outcome, method, group, label) %>%
  summarise(het_p = pchisq(sum((b - sum(b / se^2) / sum(1 / se^2))^2 / se^2), 1, lower.tail = FALSE), # Heterogeneity
            b = sum(b / se^2) / sum(1 / se^2), se = sqrt(1 / sum(1 / se^2)),
            n = sum(n), deaths = sum(deaths), .groups = "drop") %>%
  mutate(dataset = "pooled", p = 2 * pnorm(-abs(b / se)))
robust <- bind_rows(robust, pooled) %>% mutate(hr = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se))

#### 4 Forest plots ####
robust_plot <- function(m) {
  r <- robust %>%
    filter(method == m) %>%
    mutate(label = factor(label, rev(unique(sapply(specs, `[[`, "label")))),
           group = factor(group, unique(sapply(specs, `[[`, "group"))),
           outcome = factor(outcomes[outcome], outcomes),
           cohort = factor(cohort_label[dataset], cohort_label),
           sig = ifelse(p < 0.05, "p < 0.05", "p >= 0.05"))
  # position_dodge stacks the first level at the bottom; reverse so hip is on top and pooled at the bottom
  ggplot(r, aes(hr, label, colour = cohort, shape = cohort, alpha = sig, group = forcats::fct_rev(cohort))) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.6, position = position_dodge(width = 0.75)) +
    geom_point(size = 2.2, position = position_dodge(width = 0.75)) +
    facet_grid(group ~ outcome, scales = "free_y", space = "free_y") +
    scale_x_log10(breaks = c(0.5, 0.8, 1, 1.25, 1.5, 2, 3)) +
    scale_colour_manual(values = c("#2360A0", "#B74F27", "grey15")) +
    scale_shape_manual(values = c(16, 16, 18)) +
    scale_alpha_manual(values = c(1, 0.35), guide = "none") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", strip.text.y = element_text(angle = 270)) +
    labs(title = paste(if (m == "weighted") "Survey-weighted" else "Unweighted", "Cox models"),
         x = "Hazard ratio per SD of FPC2 (log scale)", y = NULL, colour = NULL, shape = NULL)
}
p_weighted <- robust_plot("weighted")
p_unweighted <- robust_plot("unweighted")
p_weighted
p_unweighted
