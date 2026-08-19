# Prospective ACL motor unit study

Reproducibility repository for the manuscript **“Quadriceps muscle sub-maximal motor unit firing characteristics and knee extensor strength following anterior cruciate ligament reconstruction.”**

## Main six-month strength analysis

The original strength-prediction analysis is reproduced by:

- `Model_DefaultPriors.R`
- `Model_StrongerPriors.R`
- `Model_WeakerPriors.R`

using:

- `ModellingData.csv`

These scripts evaluate the association between early submaximal vastii motor-unit firing characteristics and six-month knee extensor strength using Bayesian regression models.

## Secondary MUAP–firing-rate analysis

Following peer review, an additional motor-unit-level analysis was undertaken to evaluate the relationship between surface-detected motor-unit action potential (MUAP) amplitude and motor-unit firing rate, including whether this relationship differed between the ACL-reconstructed and contralateral limbs.

Files:

- `MUAP_FiringRate_Analysis.R` — complete R/brms analysis.
- `MUAP_FiringRate_Data.csv` — de-identified motor-unit-level analysis dataset.
- `MUAP_FiringRate_DataDictionary.csv` — definitions of released variables.
- `MUAP_FiringRate_ReferenceResults.csv` — archived posterior summaries from the final Stage 2F v3 crossed-hierarchical analysis used as a numerical reference.

### MUAP analysis structure

The primary analysis evaluates the change in plateau motor-unit firing rate associated with a one-standard-deviation greater log-transformed MUAP amplitude **within the same contraction**.

The analysis accounts for dependence using random intercepts for:

1. participant;
2. original Delsys motor unit within recording; and
3. contraction bout.

The secondary model includes a MUAP-by-limb interaction. Sensitivity analyses restrict observations to motor units with at least five plateau firings, progressively stricter quality-control subsets, and an alternative mean-instantaneous-firing-rate outcome.

Surface MUAP amplitude is treated as an electrophysiological characteristic and is **not** interpreted as a direct anatomical measure of motor-unit size. Recruitment threshold was not reconstructed in this analysis.

### De-identification

The public MUAP dataset uses `P###` participant identifiers. Original participant names and original recording filenames are not included. Recording identifiers were replaced with participant-specific sequential identifiers while preserving the structure required to identify repeated observations of the same decomposed motor unit across contraction bouts.

## Running the MUAP analysis

From the repository root in R:

```r
source("MUAP_FiringRate_Analysis.R")
```

Required packages:

```r
install.packages(c("brms", "tidyverse", "posterior"))
```

The script creates a `MUAP_FiringRate_outputs/` directory containing posterior summaries, limb-specific slopes, diagnostics, posterior predictive checks, and `sessionInfo()`.

The R implementation uses `brms`/Stan. The archived Stage 2F v3 reference analysis used a custom Gibbs sampler; therefore posterior values are expected to be very close but not bit-for-bit identical. The script automatically compares the R estimates against `MUAP_FiringRate_ReferenceResults.csv`.

## Data availability and DOI

This GitHub repository is archived through Zenodo. GitHub releases can therefore be used to create versioned Zenodo records/DOIs for the analysis and data accompanying the manuscript.
