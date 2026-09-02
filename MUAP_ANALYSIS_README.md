# Secondary MUAP–firing analysis

This material supports the secondary exploratory analysis of the relationship between surface-detected motor unit action potential (MUAP) amplitude and plateau motor unit firing rate reported in the revised manuscript.

## Public model-ready data

`MotorUnit_files_for_analysis/Stage2E_model_ready/Stage2E_primary_model_ready_public.csv.gz`

The public compressed CSV is a deidentified, analysis-ready version of the Stage 2E dataset. Participant and recording identifiers were replaced with stable study codes while preserving the exact grouping and repeated-measures structure used by the model. Unrelated upstream processing/QC columns were omitted; all variables required to reproduce the Stage 2F analysis were retained.

The three sensitivity datasets used in the original analysis are exact subsets of the primary Stage 2E dataset. Their eligibility flags are retained in the public CSV and the public analysis script reconstructs those subsets programmatically.

After applying the within-bout estimability requirement, the primary analysis contains 4,241 motor unit observations from 3,806 unique recording-level motor units across 510 contraction bouts and 25 participants.

## Analysis

Run from the repository root:

```bash
python stage2f_v3_reproducible_public.py
```

The analysis requires NumPy and Python's standard library. It fits:

- the primary within-bout MUAP–firing association;
- the MUAP × limb modifier model;
- sensitivity analysis requiring at least five plateau firings;
- the A-clean + B sensitivity analysis;
- the A-clean sensitivity analysis; and
- the alternative mean instantaneous firing-rate outcome.

Outputs are written to:

`MotorUnit_files_for_analysis/Stage2F_v3_final_crossed_hierarchical/`

## Figure 4

After running the Stage 2F analysis:

```bash
python figure4_muap_firing_model_based_v5.py
```

This requires NumPy and Matplotlib. The figure is written to:

`MotorUnit_files_for_analysis/Stage2F_v3_final_crossed_hierarchical/Figures/`

The solid black line is the posterior-mean model relationship. The grey ribbon represents the 95% credible interval for the limb-specific slope, anchored at the posterior-mean fitted value at standardized log-MUAP amplitude = 0.

## Interpretation

Surface MUAP amplitude is treated as an electrophysiological characteristic of the detected motor unit and is not interpreted as a direct measure of anatomical motor unit size.
