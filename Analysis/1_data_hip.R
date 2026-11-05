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
library("devtools")
if (!require("rnhanesdata")) {
  install_github("andrew-leroux/rnhanesdata")
  require("rnhanesdata")
}

#### 1 Load the data ####
# The five publicly available data categories are:
# - Demographics (DEMO)
# - Examination (EXAM)
# - Laboratory (LAB)
# - Questionnaire (Q)
# - Dietary (DIET)
# Activity count data (accelerometer)
head(PAXINTEN_C)
head(PAXINTEN_D)
# Wear/non-wear flags associated (derived from the activity count data)
head(Flags_C)
head(Flags_D)
# Mortality data linked to NHANES
mort_ls <- process_mort(mort_release_yr = 2015)
Mortality_2015_C <- mort_ls$Mortality_2015_C
Mortality_2015_D <- mort_ls$Mortality_2015_D

# NHANES demographic, survey design, lifestyle, and comorbiditiy variables
library(haven)
process_covar <- function(waves = c("C", "D"),
                          varnames = c(
                            "SDDSRVYR", "WTMEC2YR", "WTINT2YR",
                            "SDMVPSU", "SDMVSTRA",
                            "RIDAGEMN", "RIDAGEEX", "RIDAGEYR", "RIDRETH1", "RIAGENDR",
                            "BMXWT", "BMXHT", "BMXBMI", "DMDEDUC2",
                            "ALQ101", "ALQ110", "ALQ120Q", "ALQ120U", "ALQ130",
                            "SMQ020", "SMD030", "SMQ040",
                            "MCQ220", "MCQ160F", "MCQ160B", "MCQ160C",
                            "PFQ049", "PFQ054", "PFQ057", "PFQ059", "PFQ061B", "PFQ061C",
                            "DIQ010"
                            # "LBXTC", "LBXHDD", "LBDHDD",
                            ## 1. cholesterol. Note LBXHDD and LBDHDD are the same variable,
                            ## but different names for 2003-2004 and 2005-2006 cohorts
                            # "BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4"
                            ## 2. blood pressure measurements
                          ),
                          localpath = NULL,
                          extractAll = FALSE,
                          exclude_files = c()) { #### New parameter for files to exclude ####

  if (!all(waves %in% LETTERS)) stop("One or more waves invalid.")
  stopifnot(length(waves) >= 1)
  stopifnot(is.vector(waves))

  waves <- sort(waves)
  cnt <- ifelse("SEQN" %in% varnames, 1, 0)
  varnames <- unique(c("SEQN", varnames))

  if (is.null(localpath)) {
    localpath <- system.file(file.path("extdat", "covar"), package = "rnhanesdata")
    #### localpath ####
  }
  files_full <- list.files(localpath)

  # Convert exclude_files to uppercase for case-insensitive matching
  exclude_files <- toupper(exclude_files)

  if (any(!sapply(paste0("DEMO_", waves, ".XPT"), function(pattern) {
    any(grepl(pattern, files_full, ignore.case = TRUE))
  }))) { #### (case insensitive) ####
    warning(cat(paste0("One or more demographic files were not found in the data directory (", paste0("DEMO_", waves, ".XPT", collapse = ", "), ").")))
  }

  ret <- rep(list(NULL), length(waves))
  names(ret) <- paste0("Covariate_", waves)
  pb <- txtProgressBar(min = 0, max = length(waves), style = 3)
  for (i in seq_along(waves)) {
    cohort <- waves[i]

    pathExt <- paste0("_", cohort, ".", sep = "")
    # find only those files associated with a particular NHANES wave - case insensitive for XPT extension
    files <- files_full[grepl(paste0(pathExt, "xpt$", collapse = ""), files_full, ignore.case = TRUE)]

    # Exclude files that contain any of the patterns #
    files <- files[!sapply(toupper(files), function(x) any(sapply(exclude_files, function(pattern) grepl(pattern, x))))]

    if (length(files) == 0) {
      message(cat(paste0("\n No data associated with wave ", waves[i], " was found. \n")))
      setTxtProgressBar(pb, i)
      next
    }

    # find which files contain the variables requested by the user
    covarMats <- lapply(files, function(x) {
      mat <- try(read_xpt(file.path(localpath, x)))
      if (inherits(mat, "try-error")) {
        return(NULL)
      }
      if (!"SEQN" %in% colnames(mat)[1]) {
        return(NULL)
      }
      if (!extractAll) {
        if (length(setdiff(intersect(colnames(mat), varnames), "SEQN")) == 0) {
          return(NULL)
        }
        mat <- mat[, colnames(mat) %in% varnames, drop = FALSE]
      }
      if (!any(is.null(dim(mat)), all(colnames(mat) %in% "SEQN"))) {
        return(mat)
      } else {
        return(NULL)
      }
    })
    covarMats <- covarMats[!vapply(covarMats, is.null, logical(1))]
    matchedNames <- lapply(covarMats, colnames)
    numMatched <- length(unlist(matchedNames))

    if (numMatched == 0) {
      message(cat(paste0("\n No variables specified by the varnames argument was found for wave ", waves[i], "\n")))
      setTxtProgressBar(pb, i)
      next
    }
    if (numMatched > 0 & !extractAll) {
      message(
        cat(paste(
          "\n For", cohort, "cohort,",
          (numMatched - length(matchedNames) + cnt),
          "Covariates Found of", (length(varnames) - 1 + cnt), "specified.",
          "Missing covariates:",
          paste(setdiff(varnames, unlist(matchedNames)), collapse = ", ")
        ))
      )
    }

    ids <- lapply(covarMats, function(x) x[["SEQN"]])
    uids <- as.integer(sort(unique(unlist(ids))))
    rep_SEQN <- vapply(ids, function(x) any(duplicated(x)), logical(1))
    rep_inx <- which(rep_SEQN)
    notrep_inx <- which(!rep_SEQN)

    totalCols <- sum(vapply(covarMats, ncol, numeric(1))) - length(covarMats) + 1
    CovarMat <- matrix(NA, ncol = totalCols, nrow = length(uids))
    colnames(CovarMat) <- c("SEQN", unlist(sapply(covarMats, function(x) colnames(x)[-1])))
    CovarMat[, "SEQN"] <- uids
    CovarMat <- data.frame(CovarMat)
    invisible(lapply(covarMats[notrep_inx], function(x) CovarMat[, colnames(x)[-1]] <<- x[match(CovarMat$SEQN, x$SEQN), -1, drop = FALSE]))

    if (length(rep_inx) != 0) {
      invisible(lapply(covarMats[rep_inx], function(x) {
        n_vars <- length(colnames(x))
        for (j in 2:n_vars) {
          max_reps <- max(table(x$SEQN), na.rm = TRUE)
          mat_tmp <- matrix(NA, ncol = max_reps, nrow = nrow(CovarMat))

          for (i in 1:nrow(CovarMat)) {
            if (!CovarMat$SEQN[i] %in% x$SEQN) next
            ji_inx <- which(x$SEQN == CovarMat$SEQN[i])
            mat_tmp[i, 1:length(ji_inx)] <- x[ji_inx, j, drop = TRUE]

            rm(list = c("ji_inx"))
          }
          CovarMat[colnames(x)[j]] <<- I(mat_tmp)
          rm(list = c("max_reps", "mat_tmp"))
        }
      }))
      message(
        cat(paste(
          "\n Variables with repeated observations per subject found for the following variables:",
          paste(sapply(covarMats[rep_inx], function(x) colnames(x)[-1]), collapse = ","),
          "Note that these variables will be stored with class AsIs() objects in resulting data frames. ",
          "See ?I for details on AsIs class. \n"
        ))
      )
    }

    out_name <- paste0("Covariate_", cohort)
    ret[[out_name]] <- CovarMat
    rm(list = c("CovarMat", "out_name"))
    setTxtProgressBar(pb, i)
  }

  ret
}
all_covar_list <- process_covar()
Covariate_C <- all_covar_list$Covariate_C
Covariate_D <- all_covar_list$Covariate_D
# variables creation
Covariate_C$Race <- factor(Covariate_C$RIDRETH1,
  levels = 1:5,
  labels = c("Mexican American", "Other Hispanic", "White", "Black", "Other"),
  ordered = FALSE
)
Covariate_D$Race <- factor(Covariate_D$RIDRETH1,
  levels = 1:5,
  labels = c("Mexican American", "Other Hispanic", "White", "Black", "Other"),
  ordered = FALSE
)
Covariate_C$Race <- relevel(Covariate_C$Race, ref = "White")
Covariate_D$Race <- relevel(Covariate_D$Race, ref = "White")

Covariate_C$Gender <- factor(Covariate_C$RIAGENDR, levels = 1:2, labels = c("Male", "Female"), ordered = FALSE)
Covariate_D$Gender <- factor(Covariate_D$RIAGENDR, levels = 1:2, labels = c("Male", "Female"), ordered = FALSE)

Covariate_C$Diabetes <- factor(Covariate_C$DIQ010,
  levels = c(1, 2, 3, 7, 9),
  labels = c("Yes", "No", "Borderline", "Refused", "Don't know"), ordered = FALSE
)
Covariate_D$Diabetes <- factor(Covariate_D$DIQ010,
  levels = c(1, 2, 3, 7, 9),
  labels = c("Yes", "No", "Borderline", "Refused", "Don't know"), ordered = FALSE
)
Covariate_C$Diabetes <- relevel(Covariate_C$Diabetes, ref = "No")
Covariate_D$Diabetes <- relevel(Covariate_D$Diabetes, ref = "No")

Covariate_C$CHF <- factor(Covariate_C$MCQ160B,
  levels = c(1, 2, 7, 9),
  labels = c("Yes", "No", "Refused", "Don't know"), ordered = FALSE
)
Covariate_D$CHF <- factor(Covariate_D$MCQ160B,
  levels = c(1, 2, 7, 9),
  labels = c("Yes", "No", "Refused", "Don't know"), ordered = FALSE
)
Covariate_C$CHF <- relevel(Covariate_C$CHF, ref = "No")
Covariate_D$CHF <- relevel(Covariate_D$CHF, ref = "No")

Covariate_C$CHD <- factor(Covariate_C$MCQ160C,
  levels = c(1, 2, 7, 9),
  labels = c("Yes", "No", "Refused", "Don't know"), ordered = FALSE
)
Covariate_D$CHD <- factor(Covariate_D$MCQ160C,
  levels = c(1, 2, 7, 9),
  labels = c("Yes", "No", "Refused", "Don't know"), ordered = FALSE
)
Covariate_C$CHD <- relevel(Covariate_C$CHD, ref = "No")
Covariate_D$CHD <- relevel(Covariate_D$CHD, ref = "No")

Covariate_C$Cancer <- factor(Covariate_C$MCQ220,
  levels = c(1, 2, 7, 9),
  labels = c("Yes", "No", "Refused", "Don't know"), ordered = FALSE
)
Covariate_D$Cancer <- factor(Covariate_D$MCQ220,
  levels = c(1, 2, 7, 9),
  labels = c("Yes", "No", "Refused", "Don't know"), ordered = FALSE
)
Covariate_C$Cancer <- relevel(Covariate_C$Cancer, ref = "No")
Covariate_D$Cancer <- relevel(Covariate_D$Cancer, ref = "No")

Covariate_C$Stroke <- factor(Covariate_C$MCQ160F,
  levels = c(1, 2, 7, 9),
  labels = c("Yes", "No", "Refused", "Don't know"), ordered = FALSE
)
Covariate_D$Stroke <- factor(Covariate_D$MCQ160F,
  levels = c(1, 2, 7, 9),
  labels = c("Yes", "No", "Refused", "Don't know"), ordered = FALSE
)
# re-level the factor variable to have the baseline be No
Covariate_C$Stroke <- relevel(Covariate_C$Stroke, ref = "No")
Covariate_D$Stroke <- relevel(Covariate_D$Stroke, ref = "No")

Covariate_C$EducationAdult <- factor(Covariate_C$DMDEDUC2,
  levels = c(1, 2, 3, 4, 5, 7, 9),
  labels = c(
    "Less than 9th grade", "9-11th grade", "High school grad/GED or equivalent",
    "Some College or AA degree", "College graduate or above", "Refused", "Don't know"
  ),
  ordered = FALSE
)
Covariate_D$EducationAdult <- factor(Covariate_D$DMDEDUC2,
  levels = c(1, 2, 3, 4, 5, 7, 9),
  labels = c(
    "Less than 9th grade", "9-11th grade", "High school grad/GED or equivalent",
    "Some College or AA degree", "College graduate or above", "Refused", "Don't know"
  ),
  ordered = FALSE
)

Covariate_C$BMI <- Covariate_C$BMXBMI
Covariate_D$BMI <- Covariate_D$BMXBMI

Covariate_C$BMI_cat <- cut(Covariate_C$BMI,
  breaks = c(0, 18.5, 25, 30, Inf),
  labels = c("Underweight", "Normal", "Overweight", "Obese")
)
Covariate_D$BMI_cat <- cut(Covariate_D$BMI,
  breaks = c(0, 18.5, 25, 30, Inf),
  labels = c("Underweight", "Normal", "Overweight", "Obese")
)
Covariate_C$BMI_cat <- relevel(Covariate_C$BMI_cat, ref = "Normal")
Covariate_D$BMI_cat <- relevel(Covariate_D$BMI_cat, ref = "Normal")

# Creating a variable for cigarette smoking
## temporarily assign the SmokeCigs variable as individuals'
## response to SMQ040
Covariate_C$SmokeCigs <- Covariate_C$SMQ040
Covariate_D$SmokeCigs <- Covariate_D$SMQ040

## re-codes individuals who responded "No" to "Ever smoked 100 cigarettes in your life" (SMQ020)
## as -999 instead of missing (since these individuals were never asked question SMQ040)
Covariate_C$SmokeCigs[Covariate_C$SMQ020 == 2] <- -999
Covariate_D$SmokeCigs[Covariate_D$SMQ020 == 2] <- -999

## re-code individuals who answered "some days" to SMQ040 ("do you now smoke cigarettes")
## as "1". These individuals are considered current smokers by our definition
Covariate_C$SmokeCigs[Covariate_C$SmokeCigs == 2] <- 1
Covariate_D$SmokeCigs[Covariate_D$SmokeCigs == 2] <- 1

Covariate_C$SmokeCigs <- factor(Covariate_C$SmokeCigs,
  levels = c(-999, 3, 1),
  labels = c("Never", "Former", "Current"), ordered = FALSE
)
Covariate_D$SmokeCigs <- factor(Covariate_D$SmokeCigs,
  levels = c(-999, 3, 1),
  labels = c("Never", "Former", "Current"), ordered = FALSE
)

# Creating a variable for alcohol consumption
## classifies don't know/refused as missing
Covariate_C$ALQ101[Covariate_C$ALQ101 %in% c(7, 9)] <- NA
Covariate_D$ALQ101[Covariate_D$ALQ101 %in% c(7, 9)] <- NA

Covariate_C$ALQ110[Covariate_C$ALQ110 %in% c(7, 9)] <- NA
Covariate_D$ALQ110[Covariate_D$ALQ110 %in% c(7, 9)] <- NA

## get a factor variable which corresponds to "have you ever in your life had 12 drinks total over the course of a year?"
Covariate_C$Alcohol_Ever <- factor(as.numeric(Covariate_C$ALQ101 == 1 | Covariate_C$ALQ110 == 1), levels = c(1, 0), labels = c("Yes", "No"), ordered = FALSE)
Covariate_D$Alcohol_Ever <- factor(as.numeric(Covariate_D$ALQ101 == 1 | Covariate_D$ALQ110 == 1), levels = c(1, 0), labels = c("Yes", "No"), ordered = FALSE)

## re-code "how often drink alcohol over past 12 mos" = refused (777) and don't know (999) as missing
Covariate_C$ALQ120Q[Covariate_C$ALQ120Q %in% c(777, 999)] <- NA
Covariate_D$ALQ120Q[Covariate_D$ALQ120Q %in% c(777, 999)] <- NA

## re code # days drank alcohol units of refused/don't know as missing ()
## note: there are no observed values of 7/9 in these variables, but they are options in the survey
Covariate_C$ALQ120U[Covariate_C$ALQ120U %in% c(7, 9)] <- NA
Covariate_D$ALQ120U[Covariate_D$ALQ120U %in% c(7, 9)] <- NA

## re code # of drinks on those days that alcohol was drunk
## to be NA where the answer was "refused" or "dont know"
## note they changed the coding between 2003-2004 and 2005 and 2006 waves from 77/99 to 777/999
Covariate_C$ALQ130[Covariate_C$ALQ130 %in% c(77, 99)] <- NA
Covariate_D$ALQ130[Covariate_D$ALQ130 %in% c(777, 999)] <- NA

## get number of drinks per week for all individuals
multiplier <- 7 * c(1 / 7, 1 / 30, 1 / 365)

## (#days drank alcohol / unit time) * (unit time / week) * (drinks / day) = drinks/week
Covariate_C$DrinksPerWeek <- Covariate_C$ALQ120Q * (multiplier[Covariate_C$ALQ120U]) * Covariate_C$ALQ130
Covariate_D$DrinksPerWeek <- Covariate_D$ALQ120Q * (multiplier[Covariate_D$ALQ120U]) * Covariate_D$ALQ130

## recode individuals who are never drinkers OR non-drinkers  as drinking 0 drinks per week
Covariate_C$DrinksPerWeek[Covariate_C$Alcohol_Ever == "No"] <- 0
Covariate_D$DrinksPerWeek[Covariate_D$Alcohol_Ever == "No"] <- 0

Covariate_C$DrinksPerWeek[Covariate_C$ALQ120Q == 0] <- 0
Covariate_D$DrinksPerWeek[Covariate_D$ALQ120Q == 0] <- 0

## classify individuals as "non-drinker","moderate","heavy" using gender specific CDC thresholds of
##  no more than 7 drinks/week for women, and no more than  14 drinks per week for men.
cutoff <- c(14, 7)
Covariate_C$DrinkStatus <- cut(Covariate_C$DrinksPerWeek / cutoff[Covariate_C$RIAGENDR],
  breaks = c(-1, 0, 1, Inf),
  labels = c("Non-Drinker", "Moderate Drinker", "Heavy Drinker")
)
Covariate_D$DrinkStatus <- cut(Covariate_D$DrinksPerWeek / cutoff[Covariate_D$RIAGENDR],
  breaks = c(-1, 0, 1, Inf),
  labels = c("Non-Drinker", "Moderate Drinker", "Heavy Drinker")
)

## re-level the factor variable to have the baseline be moderate drinkers
Covariate_C$DrinkStatus <- relevel(Covariate_C$DrinkStatus, ref = "Moderate Drinker")
Covariate_D$DrinkStatus <- relevel(Covariate_D$DrinkStatus, ref = "Moderate Drinker")

# Creating a variable for mobility problem
Covariate_C$Difficulty_Walking <- factor(as.numeric(Covariate_C$PFQ061B == 1),
  levels = c(1, 0),
  labels = c("No Difficulty", "Any Difficulty"), ordered = FALSE
)
Covariate_D$Difficulty_Walking <- factor(as.numeric(Covariate_D$PFQ061B == 1),
  levels = c(1, 0),
  labels = c("No Difficulty", "Any Difficulty"), ordered = FALSE
)

Covariate_C$Difficulty_Stairs <- factor(as.numeric(Covariate_C$PFQ061C == 1),
  levels = c(1, 0),
  labels = c("No Difficulty", "Any Difficulty"), ordered = FALSE
)
Covariate_D$Difficulty_Stairs <- factor(as.numeric(Covariate_D$PFQ061C == 1),
  levels = c(1, 0),
  labels = c("No Difficulty", "Any Difficulty"), ordered = FALSE
)

## label anyone who requires special equipment to walk as "any difficulty"
inx_sp_equip_C <- which(Covariate_C$PFQ054 == 1)
inx_sp_equip_D <- which(Covariate_D$PFQ054 == 1)

Covariate_C$Difficulty_Walking[inx_sp_equip_C] <- "Any Difficulty"
Covariate_D$Difficulty_Walking[inx_sp_equip_D] <- "Any Difficulty"

Covariate_C$Difficulty_Stairs[inx_sp_equip_C] <- "Any Difficulty"
Covariate_D$Difficulty_Stairs[inx_sp_equip_D] <- "Any Difficulty"
rm(list = c("inx_sp_equip_C", "inx_sp_equip_D"))

# label anyone 59 and younger at interview who responds no to PFQ049, PFQ057, PFQ059 as "No difficulty"
inx_good_fn_C <- which(Covariate_C$PFQ049 == 2 & Covariate_C$PFQ057 == 2 & Covariate_C$PFQ059 == 2 & Covariate_C$RIDAGEYR <= 59)
inx_good_fn_D <- which(Covariate_D$PFQ049 == 2 & Covariate_D$PFQ057 == 2 & Covariate_D$PFQ059 == 2 & Covariate_D$RIDAGEYR <= 59)

Covariate_C$Difficulty_Walking[inx_good_fn_C] <- "No Difficulty"
Covariate_D$Difficulty_Walking[inx_good_fn_D] <- "No Difficulty"

Covariate_C$Difficulty_Stairs[inx_good_fn_C] <- "No Difficulty"
Covariate_D$Difficulty_Stairs[inx_good_fn_D] <- "No Difficulty"

Covariate_C$MobilityProblem <-
  factor(as.numeric(Covariate_C$Difficulty_Stairs == "Any Difficulty" | Covariate_C$Difficulty_Walking == "Any Difficulty"),
    levels = c(0, 1), labels = c("No Difficulty", "Any Difficulty"), ordered = FALSE
  )
Covariate_D$MobilityProblem <-
  factor(as.numeric(Covariate_D$Difficulty_Stairs == "Any Difficulty" | Covariate_D$Difficulty_Walking == "Any Difficulty"),
    levels = c(0, 1), labels = c("No Difficulty", "Any Difficulty"), ordered = FALSE
  )

#### 2 Download and process NHANES lab measurements ####
# (cholesterol, blood pressure)
# Create a (local) temporary directory
dir_tmp <- tempfile()
dir.create(dir_tmp)
if (!dir.exists(dir_tmp)) {
  dir.create(dir_tmp, showWarnings = FALSE)
}
dl_file <- function(url) {
  bn <- basename(url)
  destfile <- file.path(dir_tmp, bn)
  # if (!file.exists(destfile)) { #
  out <- download.file(url, destfile = destfile, mode = "wb")
  # } #
  stopifnot(file.exists(destfile))
}
## download the lab measurement data for the cohort 2003-2004
# Cholesterol - Total & HDL: LBXTC and LBXHDD
dl_file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/L13_C.xpt")
# Blood Pressure: BPXSY1 , BPXSY2, BPXSY3 and BPXSY4
dl_file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/BPX_C.xpt")

## download the lab measurement data for the cohort 2005-2006
# Total Cholesterol: LBXTC
dl_file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/TCHOL_D.xpt")
# HDL Cholesterol: LBDHDD
dl_file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/HDL_D.xpt")
# Blood Pressure, up to 4 measurements per person: BPXSY1 , BPXSY2, BPXSY3 and BPXSY4
dl_file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/BPX_D.xpt")

varnames <- c(
  "LBXTC", "LBXHDD", "LBDHDD",
  ## 1. cholesterol. Note LBXHDD and LBDHDD are the same variable,
  ## but different names for 2003-2004 and 2005-2006 cohorts
  "BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4"
  ## 2. blood pressure measurements
)

## load and merge the lab data
lab_data <- process_covar(varnames = varnames, localpath = dir_tmp)
## change column name for cholesterol variable that changed names
colnames(lab_data$Covariate_C)[colnames(lab_data$Covariate_C) == "LBXHDD"] <- "LBDHDD"
## combine waves
CVMarkers <- dplyr::bind_rows(lab_data$Covariate_C, lab_data$Covariate_D)
rm(list = c("lab_data", "dir_tmp", "varnames"))

#### 3 Merge the data ####
# ## re-code activity counts which are considered "non-wear" to be 0
# ## this doesn't impact many data points, most estimated non-wear times correspond to 0 counts
# PAXINTEN_C[, paste0("MIN", 1:1440)] <- PAXINTEN_C[, paste0("MIN", 1:1440)] * Flags_C[, paste0("MIN", 1:1440)]
# PAXINTEN_D[, paste0("MIN", 1:1440)] <- PAXINTEN_D[, paste0("MIN", 1:1440)] * Flags_D[, paste0("MIN", 1:1440)]

# Hip — retain wear (1); set non-wear and missing-status activity to NA
mask_hip_activity <- function(activity, flags) {
  x <- as.matrix(activity[, paste0("MIN", 1:1440), drop = FALSE])
  f <- as.matrix(flags[, paste0("MIN", 1:1440), drop = FALSE])

  x[is.na(f) | f != 1] <- NA_real_

  activity[, paste0("MIN", 1:1440)] <- x
  activity
}

PAXINTEN_C <- mask_hip_activity(PAXINTEN_C, Flags_C)
PAXINTEN_D <- mask_hip_activity(PAXINTEN_D, Flags_D)

library(dplyr)
## Merge covariate, mortality, and accelerometry data
AllAct_C <- left_join(PAXINTEN_C, Mortality_2015_C, by = "SEQN") %>%
  left_join(Covariate_C, by = c("SEQN", "SDDSRVYR"))
AllAct_D <- left_join(PAXINTEN_D, Mortality_2015_D, by = "SEQN") %>%
  left_join(Covariate_D, by = c("SEQN", "SDDSRVYR"))

AllFlags_C <- left_join(Flags_C, Mortality_2015_C, by = "SEQN") %>%
  left_join(Covariate_C, by = c("SEQN", "SDDSRVYR"))
AllFlags_D <- left_join(Flags_D, Mortality_2015_D, by = "SEQN") %>%
  left_join(Covariate_D, by = c("SEQN", "SDDSRVYR"))

## combine data for the two waves
AllAct <- bind_rows(AllAct_C, AllAct_D)
AllFlags <- bind_rows(AllFlags_C, AllFlags_D)

# merge with cardiovascular markers
AllAct <- left_join(AllAct, CVMarkers, by = "SEQN")
AllFlags <- left_join(AllFlags, CVMarkers, by = "SEQN")

## clean up the workspace again
rm(list = c("AllAct_C", "AllAct_D", "AllFlags_C", "AllFlags_D", "PAXINTEN_C", "PAXINTEN_D"))

#### 4 Covariates ####
## Code year 5 mortality,
## NAs for individuals with follow up less than 5 years and alive
AllAct$yr5_mort <- AllFlags$yr5_mort <-
  as.integer(ifelse(AllAct$permth_exm / 12 <= 5 & AllAct$mortstat == 1, 1,
    ifelse(AllAct$permth_exm / 12 < 5 & AllAct$mortstat == 0, NA, 0)
  ))

## Create Age in years using the age at examination (i.e. when participants wore the device)
AllAct$Age <- AllFlags$Age <- AllAct$RIDAGEEX / 12

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
  na.rm = T
))

# Calculate daily activity summary measures
# Decsription of these activity related measures is available at:
# https://www.biorxiv.org/content/10.1101/182337v1

## Assign just the activity and wear/non-wear flag data to matrices.
act_mat <- as.matrix(AllAct[, paste0("MIN", 1:1440)])
flag_mat <- as.matrix(AllFlags[, paste0("MIN", 1:1440)])

## replace NAs with 0s
act_mat[is.na(act_mat)] <- 0
flag_mat[is.na(flag_mat)] <- 0
act_mat[act_mat < 0] <- 0
flag_mat[flag_mat < 0] <- 0

# total activity count (TAC)
AllAct$TAC <- AllFlags$TAC <- rowSums(act_mat)
# total log activity count (TLAC)
AllAct$TLAC <- AllFlags$TLAC <- rowSums(log(1 + act_mat))
# total accelerometer wear time (WT)
AllAct$WT <- AllFlags$WT <- rowSums(flag_mat)
# total sedentary time (ST)
AllAct$ST <- AllFlags$ST <- rowSums(act_mat < 100)
# total time spent in moderate to vigorous physical activity (MVPA)
AllAct$MVPA <- AllFlags$MVPA <- rowSums(act_mat >= 2020)

## calculate fragmentation measures
bout_mat <- apply(act_mat >= 100, 1, function(x) {
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

# Calculate Peak1 (maximum single count value for each day)
AllAct$Peak1 <- AllFlags$Peak1 <- apply(act_mat, 1, function(i) {
  return(max(i, na.rm = TRUE))
})
# Calculate Peak30 (average of top 30 count values for each day)
AllAct$Peak30 <- AllFlags$Peak30 <- apply(act_mat, 1, function(i) {
  sorted_values <- sort(i, decreasing = TRUE)
  if (length(sorted_values) >= 30) {
    return(mean(sorted_values[1:30], na.rm = TRUE))
  } else {
    return(mean(sorted_values, na.rm = TRUE))
  }
})

# # compute total log activity count in each 2-hr window,
# # 2 hour (120 minutes) binning window
# tlen <- 120
# nt   <- floor(1440/tlen)
# # create a list of indices for binning into 2-hour windows
# inx_col_ls <- split(1:1440, rep(1:nt,each=tlen))
# Act_2hr    <- sapply(inx_col_ls, function(x) rowSums(log(1+act_mat[,x,drop=FALSE])))
# colnames(Act_2hr) <- paste0("TLAC_",c(1:12))
#
# AllAct   <- cbind(AllAct, Act_2hr)
# AllFlags <- cbind(AllFlags, Act_2hr)
#
# rm(list=c("tlen","nt","inx_col_ls","Act_2hr","act_mat","flag_mat","bout_mat"))

rm(list = c("act_mat", "flag_mat", "bout_mat"))

## Re-order columns so that activity and wear/non-wear flags are the last 1440 columns of data matrices.
act_cols <- which(colnames(AllAct) %in% paste0("MIN", 1:1440))
oth_cols <- which(!colnames(AllAct) %in% paste0("MIN", 1:1440))
AllAct <- AllAct[, c(oth_cols, act_cols)]
AllFlags <- AllFlags[, c(oth_cols, act_cols)]
rm(list = c("act_cols", "oth_cols"))

#### 5 Exclusion ####
table_dat <- AllAct[
  !duplicated(AllAct$SEQN),
  -which(colnames(AllAct) %in% c(
    paste0("MIN", 1:1440), "WEEKDAY",
    "TAC", "TLAC", "WT", "ST", "MVPA",
    "SBout", "ABout", "SATP", "ASTP", "Peak1", "Peak30"
  ))
]

nparticipants <- dim(table_dat)[1]

## exclusion_criteria
# exclude 85 and older at the time they wore the accelerometer (coded as NA)
# table_dat <- subset(table_dat, !(Age < 50 | is.na(Age)))
table_dat <- subset(table_dat, !is.na(Age))
nparticipants - dim(table_dat)[1]

## get the SEQN associated with individuals with fewer than *7* days accelerometer wear time
## with at least 10 hours OR had their data quality/device calibration flagged by NHANES
keep_inx <- exclude_accel(AllAct, AllFlags)
Act_Analysis <- AllAct[keep_inx, ]
Flags_Analysis <- AllFlags[keep_inx, ]
nms_rm <- unique(c(Act_Analysis$SEQN[-which(Act_Analysis$SEQN %in% names(table(Act_Analysis$SEQN))[table(Act_Analysis$SEQN) == 7])], setdiff(AllAct$SEQN, Act_Analysis$SEQN)))
# nms_rm <- unique(c(levels(Act_Analysis$SEQN)[!levels(Act_Analysis$SEQN) %in% names(which(table(Act_Analysis$SEQN) == 7))], setdiff(AllAct$SEQN, Act_Analysis$SEQN)))

rm(list = c("keep_inx"))

## Additional inclusion/exclusion criteria.
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
## Exclude = 1 indicates an individual does not meet our inclusion criteria
## Exclude = 0 indicates an individual does meet our inclusion criteria
# the expanded run
exclude_criteria <- if (expanded) criteria_vec[c(3, 4)] else criteria_vec
eval(parse(text = paste0("table_dat$Exclude <- as.integer(", paste0(exclude_criteria, collapse = "|"), ")")))

## Create our dataset for analysis with one row per subject
## containing only those subjects who meet our inclusion criteria.
data_analysis <- subset(table_dat, Exclude == 0, select = -Exclude)
max(data_analysis$permth_exm / 12)

flowchart_counts_hip <- list(
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
print(flowchart_counts_hip)  
# $n_initial
# [1] 14631

# $n_excluded_age
# [1] 300

# $n_after_age
# [1] 14331

# $n_missing_bmi
# [1] 96

# $n_missing_edu
# [1] 6001

# $n_bad_accel
# [1] 10437

# $n_missing_mort
# [1] 5029

# $n_short_followup
# [1] 0

# $n_missing_lab
# [1] 1976

# $n_step1_final
# [1] 2615

# expanded
# $n_initial
# [1] 14631

# $n_excluded_age
# [1] 300

# $n_after_age
# [1] 14331

# $n_missing_bmi
# [1] 96

# $n_missing_edu
# [1] 6001

# $n_bad_accel
# [1] 10437

# $n_missing_mort
# [1] 5029

# $n_short_followup
# [1] 0

# $n_missing_lab
# [1] 1976

# $n_step1_final
# [1] 2917
## get adjusted survey weights using the reweight_accel function
# data_analysis <- reweight_accel(data_analysis)
## Get activity/flag data for only those included participants AND days.
## Since we've already removed the "bad" days from Act_Analysis and Act_Flags,
## we need only subset based on subject ID now
Act_Analysis <- subset(Act_Analysis, SEQN %in% data_analysis$SEQN)
Flags_Analysis <- subset(Flags_Analysis, SEQN %in% data_analysis$SEQN)

# ## calculate subject specific averages of the accelerometry features
# ## using only the "good" days of data
# act_var_nms <- c("TAC","TLAC","WT","ST","MVPA","ABout","SBout","SATP","ASTP")
# for(i in act_var_nms){
#   data_analysis[[i]] <- vapply(data_analysis$SEQN, function(x) mean(Act_Analysis[[i]][Act_Analysis$SEQN==x]), numeric(1))
# }

## clean up the workspace
rm(list = c(
  "AllAct", "AllFlags", "Covariate_C", "Covariate_D",
  "criteria_vec", "nms_rm"
))

#### save ####
# Calculate percentages for count ranges: <100, 100-2020, >2020
act_mat_final <- as.matrix(Act_Analysis[, paste0("MIN", +1:1440)])
all_counts <- as.vector(act_mat_final)
count_categories <- cut(all_counts,
  breaks = c(-Inf, 100, 2020, Inf),
  labels = c("<100", "100-2020", ">2020"),
  include.lowest = TRUE
)
count_percentages <- table(count_categories) / length(all_counts) * 100

Act_Analysis$SEQN <- as.factor(Act_Analysis$SEQN)
data_analysis$SEQN <- as.factor(data_analysis$SEQN)
Flags_Analysis$SEQN <- as.factor(Flags_Analysis$SEQN)
# Convert the "weekday" column
weekdays <- c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")
Act_Analysis$WEEKDAY <- factor(Act_Analysis$WEEKDAY, levels = 1:7, labels = weekdays)
Flags_Analysis$WEEKDAY <- factor(Flags_Analysis$WEEKDAY, levels = 1:7, labels = weekdays)

Act_Analysis_subset <- Act_Analysis[, c("SEQN", "WEEKDAY", paste0("MIN", 1:1440))]
library(reshape2)
Act_Analysis_long <- melt(Act_Analysis_subset, id.vars = c("SEQN", "WEEKDAY"), variable.name = "MIN", value.name = "value")
Act_Analysis_long <- Act_Analysis_long %>%
  mutate(Accelerometry = paste(Act_Analysis_long$WEEKDAY, Act_Analysis_long$MIN, sep = "_")) %>%
  arrange(SEQN, WEEKDAY)
Act_Analysis_wide <- dcast(Act_Analysis_long, SEQN ~ factor(Accelerometry, levels = unique(Accelerometry)), value.var = "value")
rm(list = c("Act_Analysis_subset"))

# Define the directory and file paths
dir_path <- "./data/count/"
if (!dir.exists(dir_path)) {
  dir.create(dir_path)
}
# Save the object
save(Act_Analysis, file = paste0(dir_path, "Act_Analysis_old", out_suffix, ".RData"))
save(Flags_Analysis, file = paste0(dir_path, "Flags_Analysis_old", out_suffix, ".RData"))
save(Act_Analysis_long, file = paste0(dir_path, "Act_Analysis_long_old", out_suffix, ".RData"))
save(Act_Analysis_wide, file = paste0(dir_path, "Act_Analysis_wide_old", out_suffix, ".RData"))
save(data_analysis, file = paste0(dir_path, "data_analysis_old", out_suffix, ".RData"))
