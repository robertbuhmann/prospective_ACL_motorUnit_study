# Prospective ACL motor-unit study

This repository contains the de-identified data and R analysis code associated with the manuscript:

**Quadriceps muscle sub-maximal motor unit firing characteristics and knee extensor strength following anterior cruciate ligament reconstruction**

## Primary six-month strength analysis

The existing files at the repository root reproduce the Bayesian analyses examining associations between early motor-unit firing characteristics and six-month knee extensor strength:

- `ModellingData.csv`
- `Model_DefaultPriors.R`
- `Model_StrongerPriors.R`
- `Model_WeakerPriors.R`

## Secondary MUAP–firing-rate analysis

The peer-review revision added a motor-unit-level analysis examining the relationship between surface-detected motor-unit action potential (MUAP) amplitude and motor-unit firing rate.

Files:

- `MUAP_FiringRate_Data.csv` — de-identified motor-unit analytical dataset.
- `MUAP_FiringRate_DataDictionary.csv` — variable definitions.
- `MUAP_FiringRate_Analysis.R` — R/Stan implementation of the final crossed-hierarchical analysis.
- `MUAP_FiringRate_ReferenceResults.csv` — reference estimates from the final Stage 2F v3 analysis used to check the R implementation.

The primary MUAP analysis evaluates the change in firing rate associated with a 1-SD greater log-transformed MUAP amplitude **within the same contraction bout**. The model includes fixed effects for contraction intensity, muscle, limb, and timepoint, and random intercepts for participant, original recording-MU, and contraction bout. A secondary model tests whether the MUAP–firing-rate relationship differs between the ACL-reconstructed and contralateral limbs.

Sensitivity analyses restrict the data to motor units with at least five plateau firings, progressively stricter quality-control subsets, and an alternative mean-instantaneous-firing-rate outcome.

Surface MUAP amplitude is treated as an electrophysiological characteristic and is not interpreted as a direct measure of anatomical motor-unit size.

## De-identification

Participant names and original recording filenames have been removed from the public MUAP dataset. Participant identifiers use the same `P001`-style study IDs as the primary modelling dataset. Recording, bout, and motor-unit identifiers are de-identified but preserve the grouping structure required for reproducible analysis.

## Reproducibility

Run the MUAP analysis from the repository root in R:

```r
source("MUAP_FiringRate_Analysis.R")
```

The script writes model coverage, posterior summaries, limb-specific slopes, convergence diagnostics, and compact publication results to CSV files, and saves `sessionInfo()`.

The R/Stan script encodes the same final model specification used for the reported Stage 2F v3 analysis. Because independent Monte Carlo algorithms and random-number generators can produce small numerical differences, posterior estimates should agree within normal Monte Carlo error rather than match digit-for-digit.

## Software

The primary strength scripts use packages including `brms`, `mice`, `rstan`, and `tidyverse`. The MUAP analysis requires:

- R
- `rstan`
- `tidyverse`

## Data sharing

The repository contains analytical data required to reproduce the reported statistical analyses. Raw Delsys EMG exports are not included in the public repository.
