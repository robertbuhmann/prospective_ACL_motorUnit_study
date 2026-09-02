# Figure 4: model-based MUAP amplitude versus plateau firing rate by limb.
"""
Figure 4: MUAP amplitude versus plateau firing rate by limb.

Publication figure based on the final Stage 2F v3 crossed-hierarchical model.

Figure elements
---------------
- Transparent individual MU observations in the background.
- Descriptive equal-count binned means with approximate 95% CIs.
- Solid BLACK line: posterior-mean model-based MUAP/firing relationship.
- Light grey shaded band: 95% credible interval for the limb-specific slope,
  anchored at the posterior-mean fitted value at z(log MUAP)=0.
- No background grid lines.

Dependencies
------------
numpy
matplotlib

Run from the repository root:

    python "figure4_muap_firing_model_based_v5.py"

Expected inputs
---------------
MotorUnit_files_for_analysis/
  Stage2E_model_ready/
    Stage2E_primary_model_ready_public.csv.gz

  Stage2F_v3_final_crossed_hierarchical/
    Stage2F_v3_model_posterior_summary.csv
    Stage2F_v3_key_posterior_draws.csv

Outputs
-------
MotorUnit_files_for_analysis/
  Stage2F_v3_final_crossed_hierarchical/
    Figures/
      Figure4_MUAP_firing_by_limb_model_based_v5.png
      Figure4_MUAP_firing_by_limb_model_based_v5.pdf
"""
from __future__ import annotations
import csv
import gzip
import math
import statistics
from collections import defaultdict
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
ROOT = Path.cwd()
DATA_FILE = ROOT / 'MotorUnit_files_for_analysis' / 'Stage2E_model_ready' / 'Stage2E_primary_model_ready_public.csv.gz'
V3_DIR = ROOT / 'MotorUnit_files_for_analysis' / 'Stage2F_v3_final_crossed_hierarchical'
MODEL_SUMMARY_FILE = V3_DIR / 'Stage2F_v3_model_posterior_summary.csv'
POSTERIOR_DRAWS_FILE = V3_DIR / 'Stage2F_v3_key_posterior_draws.csv'
OUTPUT_DIR = V3_DIR / 'Figures'
PNG_FILE = OUTPUT_DIR / 'Figure4_MUAP_firing_by_limb_model_based_v5.png'
PDF_FILE = OUTPUT_DIR / 'Figure4_MUAP_firing_by_limb_model_based_v5.pdf'
OUTCOME_FIELD = 'plateau_train_rate_hz'
MUAP_FIELD = 'muap_peak_to_peak_mean_raw'
BOUT_FIELD = 'stage2e_bout_uid'
LIMB_FIELD = 'limb'
LIMB_MODEL_NAME = 'secondary_limb_modifier'
MUAP_SLOPE_PARAMETER = 'z_log_MUAP_within_bout'
MUAP_LIMB_INTERACTION_PARAMETER = 'zMUAP_x_LimbACL'
LIMB_MAIN_PARAMETER = 'Limb_ACL_vs_Opp'
N_BINS = 10

def safe_float(value):
    if value is None:
        return None
    text = str(value).strip()
    if text == '':
        return None
    try:
        value = float(text)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(value):
        return None
    return value

def read_csv_rows(path: Path):
    if not path.exists():
        raise FileNotFoundError(f'\nRequired file not found:\n{path}\n\nRun this script from the repository root.')
    opener = gzip.open if str(path).lower().endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def prepare_primary_model_rows(rows):
    prelim = []
    for row in rows:
        bout = str(row.get(BOUT_FIELD, '')).strip()
        limb = str(row.get(LIMB_FIELD, '')).strip()
        rate = safe_float(row.get(OUTCOME_FIELD))
        muap = safe_float(row.get(MUAP_FIELD))
        if not bout or limb not in {'Opp', 'ACL'} or rate is None or (muap is None) or (muap <= 0):
            continue
        new_row = dict(row)
        new_row['_outcome'] = rate
        new_row['_log_muap'] = math.log(muap)
        prelim.append(new_row)
    by_bout = defaultdict(list)
    for row in prelim:
        by_bout[row[BOUT_FIELD]].append(row)
    prepared = []
    for bout_rows in by_bout.values():
        log_values = [row['_log_muap'] for row in bout_rows]
        if len(log_values) < 2:
            continue
        sd_log = statistics.stdev(log_values)
        if not math.isfinite(sd_log) or sd_log <= 0:
            continue
        mean_log = statistics.mean(log_values)
        for row in bout_rows:
            row['_z_log_muap'] = (row['_log_muap'] - mean_log) / sd_log
            prepared.append(row)
    return prepared

def find_summary_row(rows, model, parameter):
    for row in rows:
        if str(row.get('model', '')).strip() == model and str(row.get('parameter', '')).strip() == parameter:
            return row
    raise RuntimeError(f"Could not find model='{model}', parameter='{parameter}' in {MODEL_SUMMARY_FILE.name}.")

def make_equal_count_bins(x, y, n_bins=10):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    order = np.argsort(x)
    x = x[order]
    y = y[order]
    groups = np.array_split(np.arange(len(x)), n_bins)
    output = []
    for idx in groups:
        if len(idx) == 0:
            continue
        bx = x[idx]
        by = y[idx]
        mean_x = float(np.mean(bx))
        mean_y = float(np.mean(by))
        if len(by) >= 2:
            se = float(np.std(by, ddof=1) / math.sqrt(len(by)))
            ci95 = 1.96 * se
        else:
            ci95 = 0.0
        output.append((mean_x, mean_y, ci95))
    return output

def load_limb_specific_slope_draws(draw_rows):
    main = {}
    interaction = {}
    for row in draw_rows:
        if str(row.get('model', '')).strip() != LIMB_MODEL_NAME:
            continue
        parameter = str(row.get('parameter', '')).strip()
        chain = str(row.get('chain', '')).strip()
        draw = str(row.get('draw', '')).strip()
        value = safe_float(row.get('value_Hz'))
        if not chain or not draw or value is None:
            continue
        key = (chain, draw)
        if parameter == MUAP_SLOPE_PARAMETER:
            main[key] = value
        elif parameter == MUAP_LIMB_INTERACTION_PARAMETER:
            interaction[key] = value
    common_keys = sorted(set(main).intersection(interaction))
    if len(common_keys) < 100:
        raise RuntimeError(f'\nCould not reconstruct sufficient paired posterior draws for the limb-specific slopes.\nMain slope draws found: {len(main)}\nInteraction draws found: {len(interaction)}\nPaired draws found: {len(common_keys)}')
    opp_draws = np.asarray([main[key] for key in common_keys], dtype=float)
    acl_draws = np.asarray([main[key] + interaction[key] for key in common_keys], dtype=float)
    return {'Opp': opp_draws, 'ACL': acl_draws}

def credible_limits(draws):
    return (float(np.mean(draws)), float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5)))

def clean_axes(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(False)
stage2e_rows = read_csv_rows(DATA_FILE)
prepared_rows = prepare_primary_model_rows(stage2e_rows)
if not prepared_rows:
    raise RuntimeError('No model-ready observations were produced.')
n_rows = len(prepared_rows)
n_bouts = len({row[BOUT_FIELD] for row in prepared_rows})
if n_rows != 4241 or n_bouts != 510:
    print(f'\nWARNING\n-----------------------------------------------\nPrepared MU observations: {n_rows} (Stage 2F v3 reported 4241)\nPrepared bouts: {n_bouts} (Stage 2F v3 reported 510)\nCheck that you are using the final Stage 2E dataset.\n')
summary_rows = read_csv_rows(MODEL_SUMMARY_FILE)
draw_rows = read_csv_rows(POSTERIOR_DRAWS_FILE)
limb_main_row = find_summary_row(summary_rows, LIMB_MODEL_NAME, LIMB_MAIN_PARAMETER)
beta_limb_mean = safe_float(limb_main_row.get('posterior_mean'))
if beta_limb_mean is None:
    raise RuntimeError('Could not read posterior_mean for Limb_ACL_vs_Opp.')
slope_draws = load_limb_specific_slope_draws(draw_rows)
opp_mean, opp_lower, opp_upper = credible_limits(slope_draws['Opp'])
acl_mean, acl_lower, acl_upper = credible_limits(slope_draws['ACL'])
overall_mean_y = float(np.mean([row['_outcome'] for row in prepared_rows]))
mean_acl_indicator = float(np.mean([1.0 if row[LIMB_FIELD] == 'ACL' else 0.0 for row in prepared_rows]))

def model_anchor(limb):
    acl_value = 1.0 if limb == 'ACL' else 0.0
    return overall_mean_y + beta_limb_mean * (acl_value - mean_acl_indicator)
panel_specs = [('Opp', 'A) Contralateral limb'), ('ACL', 'B) ACL reconstructed limb')]
fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.8), sharex=True, sharey=True)
all_x = np.asarray([row['_z_log_muap'] for row in prepared_rows], dtype=float)
all_y = np.asarray([row['_outcome'] for row in prepared_rows], dtype=float)
x_min = float(np.min(all_x))
x_max = float(np.max(all_x))
x_pad = 0.06 * (x_max - x_min)
common_xlim = (x_min - x_pad, x_max + x_pad)
y_min = min(0.0, float(np.min(all_y)))
y_max = float(np.max(all_y))
y_pad = 0.06 * (y_max - y_min)
common_ylim = (y_min - 0.25 * y_pad, y_max + y_pad)
slope_info = {'Opp': (opp_mean, opp_lower, opp_upper), 'ACL': (acl_mean, acl_lower, acl_upper)}
for ax, (limb, panel_title) in zip(axes, panel_specs):
    limb_rows = [row for row in prepared_rows if row[LIMB_FIELD] == limb]
    x = np.asarray([row['_z_log_muap'] for row in limb_rows], dtype=float)
    y = np.asarray([row['_outcome'] for row in limb_rows], dtype=float)
    points = ax.scatter(x, y, s=14, alpha=0.09, linewidths=0, rasterized=True, zorder=1)
    descriptive_colour = points.get_facecolors()[0]
    bins = make_equal_count_bins(x, y, N_BINS)
    bin_x = np.asarray([item[0] for item in bins])
    bin_y = np.asarray([item[1] for item in bins])
    bin_ci = np.asarray([item[2] for item in bins])
    ax.errorbar(bin_x, bin_y, yerr=bin_ci, fmt='o', markersize=6, linewidth=1.2, capsize=2.5, color=descriptive_colour, zorder=4)
    x_line = np.linspace(float(np.min(x)), float(np.max(x)), 300)
    anchor = model_anchor(limb)
    slope_mean, slope_lower, slope_upper = slope_info[limb]
    mean_line = anchor + slope_mean * x_line
    line_a = anchor + slope_lower * x_line
    line_b = anchor + slope_upper * x_line
    lower_line = np.minimum(line_a, line_b)
    upper_line = np.maximum(line_a, line_b)
    ax.fill_between(x_line, lower_line, upper_line, color='0.75', alpha=0.35, linewidth=0, zorder=2)
    ax.plot(x_line, mean_line, color='black', linestyle='-', linewidth=2.1, zorder=3)
    ax.set_title(panel_title, fontsize=14, pad=9)
    ax.set_xlim(*common_xlim)
    ax.set_ylim(*common_ylim)
    clean_axes(ax)
axes[0].set_ylabel('Plateau firing rate (Hz)', fontsize=12)
fig.supxlabel('log MUAP amplitude', fontsize=12, y=0.04)
fig.subplots_adjust(left=0.09, right=0.985, top=0.9, bottom=0.16, wspace=0.08)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
fig.savefig(PNG_FILE, dpi=600, bbox_inches='tight')
fig.savefig(PDF_FILE, bbox_inches='tight')
plt.close(fig)
print()
print('FIGURE 4 v5 CREATED')
print('===============================================')
print(f'MU observations plotted: {n_rows}')
print(f'Contraction bouts: {n_bouts}')
print()
print(f'Contralateral slope: {opp_mean:.4f} Hz (95% CrI {opp_lower:.4f} to {opp_upper:.4f})')
print(f'ACL reconstructed slope: {acl_mean:.4f} Hz (95% CrI {acl_lower:.4f} to {acl_upper:.4f})')
print()
print('Solid black line = posterior mean model relationship.')
print('Grey ribbon = 95% credible interval for the limb-specific slope, anchored at the posterior-mean fitted value at x=0.')
print()
print('Saved:')
print(PNG_FILE)
print(PDF_FILE)
