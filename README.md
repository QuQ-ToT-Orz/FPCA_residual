# Beyond the 24-Hour Rhythm: Daily Activity Residuals Patterns and Mortality
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23991162.svg)](https://doi.org/10.5281/zenodo.23991162)

This repository contains the analysis code and an interactive tool for a study of daily activity patterns beyond the 24-hour rest-activity rhythm. Using minute-level accelerometer data from NHANES (hip-worn activity counts, 2003-2006; wrist-worn MIMS units, 2011-2014), each participant's activity is fitted with a cosinor model, functional principal component analysis (FPCA) is applied to the cosinor residuals, and the component scores are related to cardiovascular and all-cause mortality with survey-weighted Cox models. The second residual component (FPC2), with relatively more activity in the early morning and mid-afternoon and less in the mid-morning and late evening, was associated with higher cardiovascular mortality in both cohorts.

- `Analysis/`: data preparation, cosinor and FPCA, mortality models, prediction and comparison with standard FPCA, as numbered scripts.
- `Shiny_tool/`: the interactive tool described below.

## Interactive tool

Compares one person's minute-level accelerometer data with fixed NHANES reference models.

### Run

Online: <https://quq-tot-orz.shinyapps.io/fpca-residual-tool/>

Locally, from the project root:

```r
shiny::runApp("Shiny_tool")
```

### Input

- **Activity file (`.csv` or `.rds`):** one row per day, ideally 7 days, with minute columns `MIN1` to `MIN1440` and optionally `SEQN` (ID). Use hip activity counts or wrist MIMS, to match the chosen cohort. Leave non-wear minutes empty. 
- **The demographic variables in the chosen cohort's risk model** (only age, gender and income for hip; age and income for wrist). 

### Output

- **Daily pattern:** the daily pattern score in cohort SD units and percentile, with the 24-hour rhythm's average level (MESOR), strength (amplitude) and peak time (acrophase). The score is FPC2 from the FPCA of participant-average cosinor residuals on the asinh scale.
- **Activity group:** the nearest cluster, from k-means on FPC scores of average asinh activity, with K chosen by the gap statistic. Clusters are descriptive groupings along a continuum, not distinct phenotypes.
- **Cardiovascular risk:** the 5-year risk of cardiovascular death with a 95% CI, from a survey-weighted Cox model with terms selected by cross-validated C-index.

Output is cohort-based context, not an individual diagnosis.

# Draft Manuscript

The repository contains analysis code implemented by the author using publicly available NHANES data. A research manuscript describing the methodology and findings is currently in preparation and has not yet undergone peer review.

# Funding

This project has received funding from the European Union’s Horizon Europe research and innovation programme under the Marie Skłodowska-Curie Actions grant agreement No. 101072993 (Project: LABDA).
