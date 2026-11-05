#### 0 Packages ####
gc()
rm(list = ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(ActCR)
# devtools::install_github("junruidi/ActCR")

#### 1 Load the data ####
set.seed(123)
dataset <- "mims" 
dataset <- "count"
expanded <- TRUE   # FALSE for the default run
out_suffix <- if (expanded) "_expanded" else ""
suffix <- paste0(if (dataset == "mims") "new" else "old", out_suffix)

data_path <- paste0("./data/", dataset, "/Act_Analysis_", suffix, ".RData")
load(data_path)
sample_fraction <- 0.15

sample_id <- Act_Analysis %>%
  group_by(Gender, Race, MobilityProblem) %>%
  summarise(ID = list(unique(SEQN)), .groups = "keep") %>% 
  mutate(sampled = lapply(ID, function(x) sample(x, size = round(length(x)*sample_fraction))))
sample_id <- sample_id %>%
  unnest(cols = c(sampled))

sampled_participants <- Act_Analysis %>%
  filter(SEQN %in% sample_id$sampled) %>%
  distinct(SEQN, Gender, Race, MobilityProblem) %>%
  arrange(SEQN)

Act_Analysis_subset <- Act_Analysis[
  Act_Analysis$SEQN %in% sample_id$sampled,
  c("SEQN", "WEEKDAY", paste0("MIN", 1:1440))
] %>%
  arrange(SEQN, WEEKDAY)

weekdays <- c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")
id <- c(head(sample_id, n=1)$sampled, tail(sample_id, n=1)$sampled)
id

#### 2 non-wear missing ####
# Recode masked minutes (non-wear or missing flag) to 0
Act_Analysis_mat <- as.matrix(Act_Analysis[, paste0("MIN", 1:1440)])
Act_Analysis_mat[is.na(Act_Analysis_mat)] <- 0

#### 3 asinh transformation for skewness ####
#### Inverse hyperbolic sine (IHS) transformation ####
Act_Analysis_asinh <- Act_Analysis
Act_Analysis_asinh[, paste0("MIN", 1:1440)] <- asinh(Act_Analysis_mat)

#### 4 data visualization #####
# Raw and asinh daily curves
plot_daily <- function(plot_data, weekday, data_type) {
  ggplot(plot_data[plot_data$WEEKDAY == weekday, ], aes(time, value, colour = SEQN)) +
    geom_line(linewidth = 0.5) +
    scale_x_continuous(breaks = c(1, 481, 721, 1021, 1201),
                       labels = c("12 AM", "8 AM", "12 PM", "5 PM", "8 PM")) +
    scale_colour_manual(values = c("red", "deepskyblue")) +
    theme_bw() +
    theme(legend.title = element_blank(), legend.position = "top",
          axis.text = element_text(size = 10, colour = "black", face = "bold"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(x = "", y = "", title = paste(data_type, "data: Day", weekday))
}

plot_data <- data.frame(
  SEQN = factor(rep(Act_Analysis$SEQN[which(Act_Analysis$SEQN %in% id)], each = 1440)),
  WEEKDAY = rep(Act_Analysis$WEEKDAY[which(Act_Analysis$SEQN %in% id)], each = 1440),
  time = 1:1440,
  raw = c(t(Act_Analysis_mat[which(Act_Analysis$SEQN %in% id), ])),
  asinh = c(t(as.matrix(Act_Analysis_asinh[which(Act_Analysis$SEQN %in% id), paste0("MIN", 1:1440)])))
)

plot_daily(transform(plot_data, value = raw), "Mon", "raw")
plot_daily(transform(plot_data, value = asinh), "Mon", "asinh")

#### 5 Nonparametric rhythm measures ####
# Hourly means (van Someren et al. 1999). 
# IS across the seven days of each participant; 
# IV on the concatenated hourly series.
hourly <- t(apply(Act_Analysis_mat, 1, function(x) colMeans(matrix(x, 60))))
by_subject <- lapply(split(seq_len(nrow(hourly)),
as.character(Act_Analysis$SEQN)),
function(inx) hourly[inx, , drop = FALSE])
rhythm <- data.frame(SEQN = names(by_subject),
                     IS = sapply(by_subject, IS),
                     IV = sapply(by_subject, function(x) IV(c(t(x)))),
                     row.names = NULL)

#### 6 save ####
data_dir <- paste0("./data/", dataset, "/")
save(Act_Analysis_asinh, file = paste0(data_dir, "Act_Analysis_asinh_", suffix, ".RData"))
# Add IS/IV to data_analysis
data_analysis_path <- paste0(data_dir, "data_analysis_", suffix, ".RData")
load(data_analysis_path)
data_analysis <- data_analysis %>%
  left_join(mutate(rhythm, SEQN = factor(SEQN, levels = levels(data_analysis$SEQN))), by = "SEQN")
save(data_analysis, file = data_analysis_path)
