# Beyond the 24-Hour Rhythm: Daily Activity Residuals Patterns and Mortality
[![DOI](https://zenodo.org/badge/1214625248.svg)](https://doi.org/10.5281/zenodo.23991162)
Compares one person's minute-level accelerometer data with fixed NHANES reference models.

## Run

Online: <https://quq-tot-orz.shinyapps.io/fpca-residual-tool/>

Locally, from the project root:

```r
shiny::runApp("Shiny_tool")
```

## Input

- **Activity file (`.csv` or `.rds`):** one row per day, ideally 7 days, with minute columns `MIN1` to `MIN1440` and optionally `SEQN` (ID). Use hip activity counts or wrist MIMS, to match the chosen cohort. Leave non-wear minutes empty. 
- **The demographic variables in the chosen cohort's risk model** (only age, gender and income for hip; age and income for wrist). 

## Output

- **Daily pattern:** the daily pattern score in cohort SD units and percentile, with the 24-hour rhythm's average level (MESOR), strength (amplitude) and peak time (acrophase). The score is FPC2 from the FPCA of participant-average cosinor residuals on the asinh scale.
- **Activity group:** the nearest cluster, from k-means on FPC scores of average asinh activity, with K chosen by the gap statistic. Clusters are descriptive groupings along a continuum, not distinct phenotypes.
- **Cardiovascular risk:** the probability of cardiovascular death during NHANES follow-up with a 95% CI, from the survey-weighted logistic model.

Output is cohort-based context, not an individual diagnosis.

# Funding

This project has received funding from the European Union’s Horizon Europe research and innovation programme under the Marie Skłodowska-Curie Actions grant agreement No. 101072993 (Project: LABDA).
