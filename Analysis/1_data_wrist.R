#### 0 Packages ####
gc()
rm(list = ls())
# Sensitivity run with --expanded to keep participants who lack BMI, education or lab measures 
# and survivors with less than five years of follow-up
# none of these variables enter the mortality models, and 
# excluding survivors by follow-up time depends on the outcome. 
expanded <- TRUE   # FALSE for the default run
out_suffix <- if (expanded) "_expanded" else ""

R.version # Puppy Cup
library(devtools)
# install_local("./HUO", force = TRUE, dependencies = TRUE)
library(HUO)

# Define the directory and file paths
dir_path <- "./data/mims/"
if (!dir.exists(dir_path)) {
  dir.create(dir_path)
}

#### 1 Downloading ####
# wget https://ftp.cdc.gov/pub/NHANES/LargeDataFiles/PAXMIN_H.xpt
# wget https://ftp.cdc.gov/pub/NHANES/LargeDataFiles/PAXMIN_G.xpt

# Read both files
PAXMIN_G <- read_selected_xpt("./data/mims/raw/PAXMIN_G.xpt")
PAXMIN_H <- read_selected_xpt("./data/mims/raw/PAXMIN_H.xpt")

#### 2 Load the data ####
# The five publicly available data categories are:
# - Demographics (DEMO)
# - Examination (EXAM)
# - Laboratory (LAB)
# - Questionnaire (Q)
# - Dietary (DIET)
# Activity count data (accelerometer)
library(dplyr)
library(tidyr)

G_ls <- processNHANESAccelerometer(PAXMIN_G)
H_ls <- processNHANESAccelerometer(PAXMIN_H)
PAXINTEN_G <- G_ls$accel
PAXINTEN_H <- H_ls$accel
# Wear/non-wear flags associated (derived from the activity count data)
Flags_G <- G_ls$flags
Flags_H <- H_ls$flags

load(paste(dir_path, "raw/PAXINTEN_G.RData", sep = ""))
load(paste(dir_path, "raw/PAXINTEN_H.RData", sep = ""))
load(paste(dir_path, "raw/Flags_G.RData", sep = ""))
load(paste(dir_path, "raw/Flags_H.RData", sep = ""))

# Mortality data linked to NHANES
# wget https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/NHANES_2011_2012_MORT_2019_PUBLIC.dat
# wget https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/NHANES_2013_2014_MORT_2019_PUBLIC.dat

mort_ls <- process_mort(
  waves = c("G", "H"),
  mort_release_yr = 2019
)

Mortality_2019_G <- mort_ls$Mortality_2019_G
Mortality_2019_H <- mort_ls$Mortality_2019_H
names(Mortality_2019_G)[names(Mortality_2019_G) == "seqn"] <- "SEQN"
names(Mortality_2019_H)[names(Mortality_2019_H) == "seqn"] <- "SEQN"

# Save the object
save(Mortality_2019_G, file = paste(dir_path, "Mortality_2019_G.RData", sep = ""))
save(Mortality_2019_H, file = paste(dir_path, "Mortality_2019_H.RData", sep = ""))

# NHANES demographic, survey design, lifestyle, and comorbiditiy variables
all_covar_list <- process_covar(
  extractAll = FALSE
)
# all_covar_list <- process_covar()
Covariate_G <- all_covar_list$Covariate_G
Covariate_H <- all_covar_list$Covariate_H

# variables creation
Covariate_G$Race <- factor(Covariate_G$RIDRETH1,
  levels = 1:5,
  labels = c("Mexican American", "Other Hispanic", "White", "Black", "Other"),
  ordered = FALSE
)
Covariate_H$Race <- factor(Covariate_H$RIDRETH1,
  levels = 1:5,
  labels = c("Mexican American", "Other Hispanic", "White", "Black", "Other"),
  ordered = FALSE
)
Covariate_G$Race <- relevel(Covariate_G$Race, ref = "White")
Covariate_H$Race <- relevel(Covariate_H$Race, ref = "White")

Covariate_G$Gender <- factor(Covariate_G$RIAGENDR, levels = 1:2, labels = c("Male", "Female"), ordered = FALSE)
Covariate_H$Gender <- factor(Covariate_H$RIAGENDR, levels = 1:2, labels = c("Male", "Female"), ordered = FALSE)

Covariate_G$Diabetes <- factor(Covariate_G$DIQ010,
  levels = c(1, 2, 3, 7, 9),
  labels = c("Yes", "No", "Borderline", "Refused", "Don't know"), ordered = FALSE
)
Covariate_H$Diabetes <- factor(Covariate_H$DIQ010,
  levels = c(1, 2, 3, 7, 9),
  labels = c("Yes", "No", "Borderline", "Refused", "Don't know"), ordered = FALSE
)
Covariate_G$Diabetes <- relevel(Covariate_G$Diabetes, ref = "No")
Covariate_H$Diabetes <- relevel(Covariate_H$Diabetes, ref = "No")

Covariate_G$CHF <- factor(Covariate_G$MCQ160B,
  levels = c(1, 2, 7, 9),
  labels = c("Yes", "No", "Refused", "Don't know"), ordered = FALSE
)
Covariate_H$CHF <- factor(Covariate_H$MCQ160B,
  levels = c(1, 2, 7, 9),
  labels = c("Yes", "No", "Refused", "Don't know"), ordered = FALSE
)
Covariate_G$CHF <- relevel(Covariate_G$CHF, ref = "No")
Covariate_H$CHF <- relevel(Covariate_H$CHF, ref = "No")

Covariate_G$CHD <- factor(Covariate_G$MCQ160C,
  levels = c(1, 2, 7, 9),
  labels = c("Yes", "No", "Refused", "Don't know"), ordered = FALSE
)
Covariate_H$CHD <- factor(Covariate_H$MCQ160C,
  levels = c(1, 2, 7, 9),
  labels = c("Yes", "No", "Refused", "Don't know"), ordered = FALSE
)
Covariate_G$CHD <- relevel(Covariate_G$CHD, ref = "No")
Covariate_H$CHD <- relevel(Covariate_H$CHD, ref = "No")

Covariate_G$Cancer <- factor(Covariate_G$MCQ220,
  levels = c(1, 2, 7, 9),
  labels = c("Yes", "No", "Refused", "Don't know"), ordered = FALSE
)
Covariate_H$Cancer <- factor(Covariate_H$MCQ220,
  levels = c(1, 2, 7, 9),
  labels = c("Yes", "No", "Refused", "Don't know"), ordered = FALSE
)
Covariate_G$Cancer <- relevel(Covariate_G$Cancer, ref = "No")
Covariate_H$Cancer <- relevel(Covariate_H$Cancer, ref = "No")

Covariate_G$Stroke <- factor(Covariate_G$MCQ160F,
  levels = c(1, 2, 7, 9),
  labels = c("Yes", "No", "Refused", "Don't know"), ordered = FALSE
)
Covariate_H$Stroke <- factor(Covariate_H$MCQ160F,
  levels = c(1, 2, 7, 9),
  labels = c("Yes", "No", "Refused", "Don't know"), ordered = FALSE
)
# re-level the factor variable to have the baseline be No
Covariate_G$Stroke <- relevel(Covariate_G$Stroke, ref = "No")
Covariate_H$Stroke <- relevel(Covariate_H$Stroke, ref = "No")

Covariate_G$EducationAdult <- factor(Covariate_G$DMDEDUC2,
  levels = c(1, 2, 3, 4, 5, 7, 9),
  labels = c(
    "Less than 9th grade", "9-11th grade", "High school grad/GED or equivalent",
    "Some College or AA degree", "College graduate or above", "Refused", "Don't know"
  ),
  ordered = FALSE
)
Covariate_H$EducationAdult <- factor(Covariate_H$DMDEDUC2,
  levels = c(1, 2, 3, 4, 5, 7, 9),
  labels = c(
    "Less than 9th grade", "9-11th grade", "High school grad/GED or equivalent",
    "Some College or AA degree", "College graduate or above", "Refused", "Don't know"
  ),
  ordered = FALSE
)

Covariate_G$BMI <- Covariate_G$BMXBMI
Covariate_H$BMI <- Covariate_H$BMXBMI

Covariate_G$BMI_cat <- cut(Covariate_G$BMI,
  breaks = c(0, 18.5, 25, 30, Inf),
  labels = c("Underweight", "Normal", "Overweight", "Obese")
)
Covariate_H$BMI_cat <- cut(Covariate_H$BMI,
  breaks = c(0, 18.5, 25, 30, Inf),
  labels = c("Underweight", "Normal", "Overweight", "Obese")
)
Covariate_G$BMI_cat <- relevel(Covariate_G$BMI_cat, ref = "Normal")
Covariate_H$BMI_cat <- relevel(Covariate_H$BMI_cat, ref = "Normal")

# Creating a variable for cigarette smoking
## temporarily assign the SmokeCigs variable as individuals'
## response to SMQ040
Covariate_G$SmokeCigs <- Covariate_G$SMQ040
Covariate_H$SmokeCigs <- Covariate_H$SMQ040

## re-codes individuals who responded "No" to "Ever smoked 100 cigarettes in your life" (SMQ020)
## as -999 instead of missing (since these individuals were never asked question SMQ040)
Covariate_G$SmokeCigs[Covariate_G$SMQ020 == 2] <- -999
Covariate_H$SmokeCigs[Covariate_H$SMQ020 == 2] <- -999

## re-code individuals who answered "some days" to SMQ040 ("do you now smoke cigarettes")
## as "1". These individuals are considered current smokers by our definition
Covariate_G$SmokeCigs[Covariate_G$SmokeCigs == 2] <- 1
Covariate_H$SmokeCigs[Covariate_H$SmokeCigs == 2] <- 1

## finally, create the factor variable based on our re-coding
Covariate_G$SmokeCigs <- factor(Covariate_G$SmokeCigs,
  levels = c(-999, 3, 1),
  labels = c("Never", "Former", "Current"), ordered = FALSE
)
Covariate_H$SmokeCigs <- factor(Covariate_H$SmokeCigs,
  levels = c(-999, 3, 1),
  labels = c("Never", "Former", "Current"), ordered = FALSE
)

# Creating a variable for alcohol consumption
## classifies don't know/refused as missing
Covariate_G$ALQ101[Covariate_G$ALQ101 %in% c(7, 9)] <- NA
Covariate_H$ALQ101[Covariate_H$ALQ101 %in% c(7, 9)] <- NA

Covariate_G$ALQ110[Covariate_G$ALQ110 %in% c(7, 9)] <- NA
Covariate_H$ALQ110[Covariate_H$ALQ110 %in% c(7, 9)] <- NA

## get a factor variable which corresponds to "have you ever in your life had 12 drinks total over the course of a year?"
Covariate_G$Alcohol_Ever <- factor(as.numeric(Covariate_G$ALQ101 == 1 | Covariate_G$ALQ110 == 1), levels = c(1, 0), labels = c("Yes", "No"), ordered = FALSE)
Covariate_H$Alcohol_Ever <- factor(as.numeric(Covariate_H$ALQ101 == 1 | Covariate_H$ALQ110 == 1), levels = c(1, 0), labels = c("Yes", "No"), ordered = FALSE)

## re-code "how often drink alcohol over past 12 mos" = refused (777) and don't know (999) as missing
Covariate_G$ALQ120Q[Covariate_G$ALQ120Q %in% c(777, 999)] <- NA
Covariate_H$ALQ120Q[Covariate_H$ALQ120Q %in% c(777, 999)] <- NA

## re code # days drank alcohol units of refused/don't know as missing ()
## note: there are no observed values of 7/9 in these variables, but they are options in the survey
Covariate_G$ALQ120U[Covariate_G$ALQ120U %in% c(7, 9)] <- NA
Covariate_H$ALQ120U[Covariate_H$ALQ120U %in% c(7, 9)] <- NA

## re code # of drinks on those days that alcohol was drunk
## to be NA where the answer was "refused" or "dont know"
## note they changed the coding between 2003-2004 and 2005 and 2006 waves from 77/99 to 777/999
Covariate_G$ALQ130[Covariate_G$ALQ130 %in% c(77, 99)] <- NA
Covariate_H$ALQ130[Covariate_H$ALQ130 %in% c(777, 999)] <- NA

## get number of drinks per week for all individuals
multiplier <- 7 * c(1 / 7, 1 / 30, 1 / 365)

## (#days drank alcohol / unit time) * (unit time / week) * (drinks / day) = drinks/week
Covariate_G$DrinksPerWeek <- Covariate_G$ALQ120Q * (multiplier[Covariate_G$ALQ120U]) * Covariate_G$ALQ130
Covariate_H$DrinksPerWeek <- Covariate_H$ALQ120Q * (multiplier[Covariate_H$ALQ120U]) * Covariate_H$ALQ130

## recode individuals who are never drinkers OR non-drinkers  as drinking 0 drinks per week
Covariate_G$DrinksPerWeek[Covariate_G$Alcohol_Ever == "No"] <- 0
Covariate_H$DrinksPerWeek[Covariate_H$Alcohol_Ever == "No"] <- 0

Covariate_G$DrinksPerWeek[Covariate_G$ALQ120Q == 0] <- 0
Covariate_H$DrinksPerWeek[Covariate_H$ALQ120Q == 0] <- 0

## classify individuals as "non-drinker","moderate","heavy" using gender specific CDC thresholds of
##  no more than 7 drinks/week for women, and no more than  14 drinks per week for men.
##  note we do not have
cutoff <- c(14, 7)
Covariate_G$DrinkStatus <- cut(Covariate_G$DrinksPerWeek / cutoff[Covariate_G$RIAGENDR],
  breaks = c(-1, 0, 1, Inf),
  labels = c("Non-Drinker", "Moderate Drinker", "Heavy Drinker")
)
Covariate_H$DrinkStatus <- cut(Covariate_H$DrinksPerWeek / cutoff[Covariate_H$RIAGENDR],
  breaks = c(-1, 0, 1, Inf),
  labels = c("Non-Drinker", "Moderate Drinker", "Heavy Drinker")
)

## re-level the factor variable to have the baseline be moderate drinkers
Covariate_G$DrinkStatus <- relevel(Covariate_G$DrinkStatus, ref = "Moderate Drinker")
Covariate_H$DrinkStatus <- relevel(Covariate_H$DrinkStatus, ref = "Moderate Drinker")

# Creating a variable for mobility problem
Covariate_G$Difficulty_Walking <- factor(as.numeric(Covariate_G$PFQ061B == 1),
  levels = c(1, 0),
  labels = c("No Difficulty", "Any Difficulty"), ordered = FALSE
)
Covariate_H$Difficulty_Walking <- factor(as.numeric(Covariate_H$PFQ061B == 1),
  levels = c(1, 0),
  labels = c("No Difficulty", "Any Difficulty"), ordered = FALSE
)

Covariate_G$Difficulty_Stairs <- factor(as.numeric(Covariate_G$PFQ061C == 1),
  levels = c(1, 0),
  labels = c("No Difficulty", "Any Difficulty"), ordered = FALSE
)
Covariate_H$Difficulty_Stairs <- factor(as.numeric(Covariate_H$PFQ061C == 1),
  levels = c(1, 0),
  labels = c("No Difficulty", "Any Difficulty"), ordered = FALSE
)

## label anyone who requires special equipment to walk as "any difficulty"
inx_sp_equip_G <- which(Covariate_G$PFQ054 == 1)
inx_sp_equip_H <- which(Covariate_H$PFQ054 == 1)

Covariate_G$Difficulty_Walking[inx_sp_equip_G] <- "Any Difficulty"
Covariate_H$Difficulty_Walking[inx_sp_equip_H] <- "Any Difficulty"

Covariate_G$Difficulty_Stairs[inx_sp_equip_G] <- "Any Difficulty"
Covariate_H$Difficulty_Stairs[inx_sp_equip_H] <- "Any Difficulty"
rm(list = c("inx_sp_equip_G", "inx_sp_equip_H"))

# label anyone 59 and younger at interview who responds no to PFQ049, PFQ057, PFQ059 as "No difficulty"
inx_good_fn_G <- which(Covariate_G$PFQ049 == 2 & Covariate_G$PFQ057 == 2 & Covariate_G$PFQ059 == 2 & Covariate_G$RIDAGEYR <= 59)
inx_good_fn_H <- which(Covariate_H$PFQ049 == 2 & Covariate_H$PFQ057 == 2 & Covariate_H$PFQ059 == 2 & Covariate_H$RIDAGEYR <= 59)

Covariate_G$Difficulty_Walking[inx_good_fn_G] <- "No Difficulty"
Covariate_H$Difficulty_Walking[inx_good_fn_H] <- "No Difficulty"

Covariate_G$Difficulty_Stairs[inx_good_fn_G] <- "No Difficulty"
Covariate_H$Difficulty_Stairs[inx_good_fn_H] <- "No Difficulty"

Covariate_G$MobilityProblem <-
  factor(as.numeric(Covariate_G$Difficulty_Stairs == "Any Difficulty" | Covariate_G$Difficulty_Walking == "Any Difficulty"),
    levels = c(0, 1), labels = c("No Difficulty", "Any Difficulty"), ordered = FALSE
  )
Covariate_H$MobilityProblem <-
  factor(as.numeric(Covariate_H$Difficulty_Stairs == "Any Difficulty" | Covariate_H$Difficulty_Walking == "Any Difficulty"),
    levels = c(0, 1), labels = c("No Difficulty", "Any Difficulty"), ordered = FALSE
  )

save(Covariate_G, file = paste(dir_path, "Covariate_G.RData", sep = ""))
save(Covariate_H, file = paste(dir_path, "Covariate_H.RData", sep = ""))

#### 3 Download and process NHANES lab measurements ####
lab2011 <- extract_nhanes_data(2011, exclude_tables = "PAHS_G", extractAll = FALSE)
lab2013 <- extract_nhanes_data(2013, exclude_tables = "PAHS_H", extractAll = FALSE)

colnames(lab2011$data)[colnames(lab2011$data) == "TCHOL_G.LBXTC"] <- "LBXTC"
colnames(lab2011$data)[colnames(lab2011$data) == "HDL_G.LBDHDD"] <- "LBDHDD"
colnames(lab2013$data)[colnames(lab2013$data) == "TCHOL_H.LBXTC"] <- "LBXTC"
colnames(lab2013$data)[colnames(lab2013$data) == "HDL_H.LBDHDD"] <- "LBDHDD"

save(lab2011, file = paste(dir_path, "lab2011.RData", sep = ""))
save(lab2013, file = paste(dir_path, "lab2013.RData", sep = ""))

load(paste(dir_path, "lab2011.RData", sep = ""))
load(paste(dir_path, "lab2013.RData", sep = ""))
CVMarkers <- dplyr::bind_rows(lab2011$data, lab2013$data)

# exclude_missing <- function(df, threshold = 0.5, essential_vars = NULL) {
#   # Calculate missing summary
#   missing_summary <- data.frame(
#     variable = names(df),
#     missing_count = colSums(is.na(df)),
#     missing_prop = colMeans(is.na(df))
#   )
#   missing_summary <- missing_summary[order(-missing_summary$missing_prop), ]
#   # Handle essential variables if specified
#   if (!is.null(essential_vars)) {
#     vars_to_check <- setdiff(names(df), essential_vars)
#     keep_vars <- c(
#       intersect(essential_vars, names(df)),
#       vars_to_check[colMeans(is.na(df[vars_to_check])) <= threshold]
#     )
#     return(list(
#       data = df[, keep_vars],
#       summary = missing_summary
#     ))
#   }
#   # Simple threshold-based exclusion if no essential variables
#   return(list(
#     data = df[, colMeans(is.na(df)) <= threshold],
#     summary = missing_summary
#   ))
# }

rm(list = c("inx_good_fn_G", "inx_good_fn_H", "lab2011", "lab2013"))

#### 4 Merge the data ####
# ## re-code activity counts which are considered "non-wear" and "sleep wear" to be 0
# # Create binary flags (1 for valid wake wear time, 0 for non-wear)
# # Convert 2/3/4 to 0 and keep 1 as 1
# # This was different in the second paper
# PAXINTEN_G[, paste0("MIN", 1:1440)] <- PAXINTEN_G[, paste0("MIN", 1:1440)] * (Flags_G[, paste0("MIN", 1:1440)] == 1)
# PAXINTEN_H[, paste0("MIN", 1:1440)] <- PAXINTEN_H[, paste0("MIN", 1:1440)] * (Flags_H[, paste0("MIN", 1:1440)] == 1)

# Wrist — retain wake wear (1) and sleep wear (2); set non-wear (3), unknown (4), and missing-status activity to NA

mask_wrist_activity <- function(activity, flags) {
  x <- as.matrix(activity[, paste0("MIN", 1:1440), drop = FALSE])
  f <- as.matrix(flags[, paste0("MIN", 1:1440), drop = FALSE])

  x[is.na(f) | !(f == 1 | f == 2)] <- NA_real_

  activity[, paste0("MIN", 1:1440)] <- x
  activity
}

PAXINTEN_G <- mask_wrist_activity(PAXINTEN_G, Flags_G)
PAXINTEN_H <- mask_wrist_activity(PAXINTEN_H, Flags_H)

library(dplyr)
## Merge covariate, mortality, and accelerometry data
AllAct_G <- left_join(PAXINTEN_G, Mortality_2019_G, by = "SEQN") %>%
  left_join(Covariate_G, by = "SEQN")
 # Remove the duplicate 
Mortality_2019_H <- Mortality_2019_H %>% distinct(SEQN, .keep_all = T)
AllAct_H <- left_join(PAXINTEN_H, Mortality_2019_H, by = "SEQN") %>%
  left_join(Covariate_H, by = "SEQN")

AllFlags_G <- left_join(Flags_G, Mortality_2019_G, by = "SEQN") %>%
  left_join(Covariate_G, by = "SEQN")
AllFlags_H <- left_join(Flags_H, Mortality_2019_H, by = "SEQN") %>%
  left_join(Covariate_H, by = "SEQN", relationship = "many-to-many")

## combine data for the two waves
AllAct <- bind_rows(AllAct_G, AllAct_H)
AllFlags <- bind_rows(AllFlags_G, AllFlags_H)

# merge with cardiovascular markers
AllAct <- left_join(AllAct, CVMarkers, by = "SEQN")
AllFlags <- left_join(AllFlags, CVMarkers, by = "SEQN")

## clean up the workspace again
rm(list = c("AllAct_G", "AllAct_H", "AllFlags_G", "AllFlags_H", "PAXINTEN_G", "PAXINTEN_H"))

#### 5 Covariates ####
## Code year 5 mortality,
## NAs for individuals with follow up less than 5 years and alive
AllAct$yr5_mort <- AllFlags$yr5_mort <-
  as.integer(ifelse(AllAct$permth_exm / 12 <= 5 & AllAct$mortstat == 1, 1,
    ifelse(AllAct$permth_exm / 12 < 5 & AllAct$mortstat == 0, NA, 0)
  ))

## Create Age in years using the age at examination (i.e. when participants wore the device)
# AllAct$Age <- AllFlags$Age <- AllAct$RIDAGEEX / 12
AllAct$Age <- AllFlags$Age <- AllAct$RIDAGEYR

## Re-level comorbidities to assign refused/don't know as not having the condition
## Note that in practice this does not affect many individuals, but it is an assumption we're making.
levels(AllAct$CHD) <- levels(AllFlags$CHD) <- list("No" = c("No", "Refused", "Don't know"), "Yes" = c("Yes"))
levels(AllAct$CHF) <- levels(AllFlags$CHF) <- list("No" = c("No", "Refused", "Don't know"), "Yes" = c("Yes"))
levels(AllAct$Stroke) <- levels(AllFlags$Stroke) <- list("No" = c("No", "Refused", "Don't know"), "Yes" = c("Yes"))
levels(AllAct$Cancer) <- levels(AllFlags$Cancer) <- list("No" = c("No", "Refused", "Don't know"), "Yes" = c("Yes"))
levels(AllAct$Diabetes) <- levels(AllFlags$Diabetes) <- list("No" = c("No", "Borderline", "Refused", "Don't know"), "Yes" = c("Yes"))

## Combine race/ethnicity: Mexican American + Other Hispanic = Hispanic
levels(AllAct$Race) <- levels(AllFlags$Race) <- list(
  "White" = "White",
  "Hispanic" = c("Mexican American", "Other Hispanic"),
  "Black" = "Black",
  "Other" = "Other"
)

## Re-level education to 2 levels (high school or less vs more) and categorize don't know/refused to be missing
levels(AllAct$EducationAdult) <- levels(AllFlags$EducationAdult) <- list(
  "High school or less" = c("Less than 9th grade", "9-11th grade", "High school grad/GED or equivalent"),
  "More than high school" = c("Some College or AA degree", "College graduate or above")
)

## Re-level alcohol consumption to include a level for "missing"
levels(AllAct$DrinkStatus) <- levels(AllFlags$DrinkStatus) <- c(
  levels(AllAct$DrinkStatus),
  "Missing"
)
AllAct$DrinkStatus[is.na(AllAct$DrinkStatus)] <-
  AllFlags$DrinkStatus[is.na(AllAct$DrinkStatus)] <- "Missing"

## Collapse alcohol to 2 levels: current drinker (moderate or heavy) vs not (non-drinker or Missing)
levels(AllAct$DrinkStatus) <- levels(AllFlags$DrinkStatus) <- list(
  "Not current drinker" = c("Non-Drinker", "Missing"),
  "Current drinker" = c("Moderate Drinker", "Heavy Drinker")
)

# average systolic blood pressure calculation
AllAct$SYS <- AllFlags$SYS <- round(apply(AllAct[, c("BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4")],
  1, mean,
  na.rm = TRUE
))

# Calculate daily activity summary measures
# Decsription of these activity related measures is available at:
# https://doi-org.vu-nl.idm.oclc.org/10.1111/sms.14762

## Assign just the activity and wear/non-wear flag data to matrices.
act_mat <- as.matrix(AllAct[, paste0("MIN", 1:1440)])
flag_mat <- as.matrix(AllFlags[, paste0("MIN", 1:1440)])

## replace NAs with 0s
act_mat[is.na(act_mat)] <- 0
flag_mat[is.na(flag_mat)] <- 0
act_mat[act_mat < 0] <- 0
flag_mat[flag_mat < 0] <- 0

# Calculate MIMS cutoffs based on percentages: 75%, 23.5%, 1.5%
cutoff_75 <- quantile(act_mat, 0.75)
cutoff_98.5 <- quantile(act_mat, 0.985)
sedentary_cutoff <- 10.558
mvpa_cutoff <- 19.614

# total activity MIMS (TAC)
AllAct$TAC <- AllFlags$TAC <- rowSums(act_mat)
# total accelerometer sleep wear time (SPT)
AllAct$SPT <- AllFlags$SPT <- rowSums(flag_mat == 2, na.rm = TRUE)
# total accelerometer wake wear time (WT)
AllAct$WT <- AllFlags$WT <- rowSums(flag_mat == 1, na.rm = TRUE)
# total sedentary time (ST)
AllAct$ST <- AllFlags$ST <- rowSums(act_mat < sedentary_cutoff)
# total time spent in moderate to vigorous physical activity (MVPA)
AllAct$MVPA <- AllFlags$MVPA <- rowSums(act_mat >= mvpa_cutoff)

## calculate fragmentation measures
bout_mat <- apply(act_mat >= sedentary_cutoff, 1, function(x) {
  mat <- rle(x) # Run Length
  sed <- mat$lengths[which(mat$values == FALSE)] # the lengths of runs of 'sedentary'
  act <- mat$lengths[mat$values == TRUE] # the lengths of runs of 'active'

  # average bout duration
  sed <- ifelse(length(sed) == 0, NA, mean(sed))
  act <- ifelse(length(act) == 0, NA, mean(act))
  c(sed, act)
})

# average bout duratio
AllAct$SBout <- AllFlags$SBout <- bout_mat[1, ]
AllAct$ABout <- AllFlags$ABout <- bout_mat[2, ]
# the reciprocal of average bout duration
AllAct$SATP <- AllFlags$SATP <- 1 / AllAct$SBout
AllAct$ASTP <- AllFlags$ASTP <- 1 / AllAct$ABout

# Calculate Peak1 (maximum single MIMS value for each day)
AllAct$Peak1 <- AllFlags$Peak1 <- apply(act_mat, 1, function(i) {
  return(max(i, na.rm = TRUE))
})
# Calculate Peak30 (average of top 30 MIMS values for each day)
AllAct$Peak30 <- AllFlags$Peak30 <- apply(act_mat, 1, function(i) {
  sorted_values <- sort(i, decreasing = TRUE)
  if (length(sorted_values) >= 30) {
    return(mean(sorted_values[1:30], na.rm = TRUE))
  } else {
    return(mean(sorted_values, na.rm = TRUE))
  }
})

rm(list = c("act_mat", "flag_mat"))

## Re-order columns so that activity and wear/non-wear flags are the last 1440 columns of data matrices.
act_cols <- which(colnames(AllAct) %in% paste0("MIN", 1:1440))
oth_cols <- which(!colnames(AllAct) %in% paste0("MIN", 1:1440))
AllAct <- AllAct[, c(oth_cols, act_cols)]
AllFlags <- AllFlags[, c(oth_cols, act_cols)]
rm(list = c("act_cols", "oth_cols"))

#### 6 Exclusion ####
table_dat <- AllAct[
  !duplicated(AllAct$SEQN),
  -which(colnames(AllAct) %in% c(
    paste0("MIN", 1:1440), "WEEKDAY",
    "TAC", "SPT", "WT", "ST", "MVPA",
    "SBout", "ABout", "SATP", "ASTP", "Peak1", "Peak30"
  ))
]

nparticipants <- dim(table_dat)[1] # a total of 6905 participants

## exclusion_criteria
# exclude 85 and older at the time they wore the accelerometer (coded as NA)
# table_dat <- subset(table_dat, !(Age < 50 | is.na(Age)))
table_dat <- subset(table_dat, !is.na(Age))
nparticipants - dim(table_dat)[1]

## get the SEQN associated with individuals with fewer than *7* days accelerometer wear time
## with at least 10 hours OR had their data quality/device calibration flagged by NHANES
keep_inx <- check_valid_wake_days(AllAct, AllFlags)
Act_Analysis <- AllAct[keep_inx, ]
Flags_Analysis <- AllFlags[keep_inx, ]
nms_rm <- unique(c(
  Act_Analysis$SEQN[-which(Act_Analysis$SEQN %in%
    names(table(Act_Analysis$SEQN))[table(Act_Analysis$SEQN) == 7])],
  setdiff(AllAct$SEQN, Act_Analysis$SEQN)
))
# nms_rm <- unique(c(
#   levels(Act_Analysis$SEQN)[!levels(Act_Analysis$SEQN) %in% names(which(table(Act_Analysis$SEQN) == 7))],
#   setdiff(AllAct$SEQN, Act_Analysis$SEQN)
# ))

rm(list = c("keep_inx"))

## Additional inclusion/exclusion criteria
criteria_vec <- c(
  "(is.na(table_dat$BMI_cat))", # missing BMI
  "(is.na(table_dat$EducationAdult))", # missing education
  "(table_dat$SEQN %in% nms_rm)", # too few "good" days of accelerometery data
  "((!table_dat$eligstat %in% 1) | is.na(table_dat$mortstat) | is.na(table_dat$permth_exm))",
  # missing mortality data (accidental deaths are kept)
  "(table_dat$mortstat == 0 & table_dat$permth_exm/12 < 5)",
  # less than 5 years of follow up with no mortality

  "(is.na(table_dat$SYS) | (is.na(table_dat$LBXTC)) | (is.na(table_dat$LBDHDD)) )"
  # missing lab measures
)

## create matrix of pairwise missing data based on our exclusion criteria
tab_miss <- matrix(NA, ncol = length(criteria_vec), nrow = length(criteria_vec))
for (i in seq_along(criteria_vec)) {
  for (j in seq_along(criteria_vec)) {
    eval(parse(text = paste0("miss_cur <- which(", criteria_vec[i], "&", criteria_vec[j], ")")))
    tab_miss[i, j] <- length(miss_cur)
    rm(list = c("miss_cur"))
  }
}
names_spaced <- c("BMI", "Education", "Bad Accel Data", "Mortality", "Follow-up", "Lab")
rownames(tab_miss) <- colnames(tab_miss) <- names_spaced

# missing BMI or education predictor variables :
tab_miss[1, 1] + tab_miss[2, 2]
# had fewer than 7 days of data with at least 10 hours of estimated wear time
tab_miss[3, 3]
# missing mortality information
tab_miss[4, 4]
# alive with follow up less than 5 years
tab_miss[5, 5]
# missing systolic blood pressure, total cholesterol or HDL cholesterol  measurements
tab_miss[6, 6]

## add in column indicating exclusion:
##   Exclude = 1 indicates an individual does not meet our inclusion criteria
##   Exclude = 0 indicates an individual does meet our inclusion criteria
# the expanded run
exclude_criteria <- if (expanded) criteria_vec[c(3, 4)] else criteria_vec
eval(parse(text = paste0("table_dat$Exclude <- as.integer(", paste0(exclude_criteria, collapse = "|"), ")")))

## Create our dataset for analysis with one row per subject
## containing only those subjects who meet our inclusion criteria.
data_analysis <- subset(table_dat, Exclude == 0, select = -Exclude)
max(data_analysis$permth_exm / 12)

flowchart_counts_wrist <- list(
  n_initial = nparticipants,
  n_excluded_age = nparticipants - nrow(table_dat),
  n_after_age = nrow(table_dat),
  n_missing_bmi = sum(is.na(table_dat$BMI_cat)),
  n_missing_edu = sum(is.na(table_dat$EducationAdult)),
  n_bad_accel = sum(table_dat$SEQN %in% nms_rm),
  n_missing_mort = sum(
    (!table_dat$eligstat %in% 1) |
      is.na(table_dat$mortstat) |
      is.na(table_dat$permth_exm),
    na.rm = TRUE
  ),
  n_short_followup = sum(
    table_dat$mortstat == 0 & (table_dat$permth_exm / 12 < 5),
    na.rm = TRUE
  ),
  n_missing_lab = sum(
    is.na(table_dat$SYS) | is.na(table_dat$LBXTC) | is.na(table_dat$LBDHDD),
    na.rm = TRUE
  ),
  n_step1_final = nrow(data_analysis)
)
print(flowchart_counts_wrist)
# $n_initial
# [1] 6905

# $n_excluded_age
# [1] 0

# $n_after_age
# [1] 6905

# $n_missing_bmi
# [1] 82

# $n_missing_edu
# [1] 1564

# $n_bad_accel
# [1] 2554

# $n_missing_mort
# [1] 1318

# $n_short_followup
# [1] 18

# $n_missing_lab
# [1] 932

# $n_step1_final
# [1] 3208

# expanded
# $n_initial
# [1] 6905

# $n_excluded_age
# [1] 0

# $n_after_age
# [1] 6905

# $n_missing_bmi
# [1] 82

# $n_missing_edu
# [1] 1564

# $n_bad_accel
# [1] 2554

# $n_missing_mort
# [1] 1318

# $n_short_followup
# [1] 18

# $n_missing_lab
# [1] 932

# $n_step1_final
# [1] 3605

## get adjusted survey weights using the reweight_accel function
# data_analysis  <- reweight_accel(data_analysis)
## Get activity/flag data for only those included participants AND days.
## Since we've already removed the "bad" days from Act_Analysis and Act_Flags,
## we need only subset based on subject ID now
Act_Analysis <- subset(Act_Analysis, SEQN %in% data_analysis$SEQN)
Flags_Analysis <- subset(Flags_Analysis, SEQN %in% data_analysis$SEQN)

## clean up the workspace
rm(list = c(
  "AllAct", "AllFlags", "Covariate_G", "Covariate_H",
  "criteria_vec", "nms_rm"
))

#### save ####
# Calculate percentages for count ranges: <100, 100-2020, >2020
act_mat_final <- as.matrix(Act_Analysis[, paste0("MIN", 1:1440)])
all_counts <- as.vector(act_mat_final)
count_categories <- cut(all_counts,
  breaks = c(-Inf, sedentary_cutoff, mvpa_cutoff, Inf),
  include.lowest = TRUE
)
count_percentages <- table(count_categories) / length(all_counts) * 100

Act_Analysis$SEQN <- as.factor(Act_Analysis$SEQN)
data_analysis$SEQN <- as.factor(data_analysis$SEQN)
Flags_Analysis$SEQN <- as.factor(Flags_Analysis$SEQN)
# Convert the "weekday" column
weekdays <- c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")
Act_Analysis$WEEKDAY <- factor(Act_Analysis$PAXDAYWM, levels = 1:7, labels = weekdays)
Flags_Analysis$WEEKDAY <- factor(Flags_Analysis$PAXDAYWM, levels = 1:7, labels = weekdays)

Act_Analysis_subset <- Act_Analysis[, c("SEQN", "WEEKDAY", paste0("MIN", 1:1440))]
library(reshape2)
Act_Analysis_long <- melt(Act_Analysis_subset, id.vars = c("SEQN", "WEEKDAY"), variable.name = "MIN", value.name = "value")
Act_Analysis_long <- Act_Analysis_long %>%
  mutate(Accelerometry = paste(Act_Analysis_long$WEEKDAY, Act_Analysis_long$MIN, sep = "_")) %>%
  arrange(SEQN, WEEKDAY)
Act_Analysis_wide <- reshape2::dcast(Act_Analysis_long, SEQN ~ factor(Accelerometry, levels = unique(Accelerometry)), value.var = "value")
rm(list = c("Act_Analysis_subset"))

save(Act_Analysis, file = paste0(dir_path, "Act_Analysis_new", out_suffix, ".RData"))
save(Flags_Analysis, file = paste0(dir_path, "Flags_Analysis_new", out_suffix, ".RData"))
save(Act_Analysis_long, file = paste0(dir_path, "Act_Analysis_long_new", out_suffix, ".RData"))
save(Act_Analysis_wide, file = paste0(dir_path, "Act_Analysis_wide_new", out_suffix, ".RData"))
save(data_analysis, file = paste0(dir_path, "data_analysis_new", out_suffix, ".RData"))
