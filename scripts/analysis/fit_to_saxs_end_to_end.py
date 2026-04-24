import numpy as np
import matplotlib.pyplot as plt
import pickle
import tqdm
import os
from scipy.spatial.distance import pdist
import argparse

parser = argparse.ArgumentParser(
    description="SAXS end to end Dmax Validation"
)
parser.add_argument("--rmf", type=str, required=True, help="Path to the combined RMF file.")
parser.add_argument("--output_dir", type=str, required=True, help="Path to the parent validation directory.")
args = parser.parse_args()

combined_rmf = args.rmf
save_path    = os.path.join(args.output_dir, "saxs_end_to_end_validation")
os.makedirs(save_path, exist_ok=True)

# Your specific residue definitions
saxs_dmax = [
    {"prot": "cipa", "res_start": 1,    "res_end": 1853, "target_distance": 575.0,  "threshold": 20.0, "label": "cellulosome end to end"}
]

# Load data
with open('./extracted_xyzr/saved_data', 'rb') as f:
    m, p = pickle.load(f)
with open('./contact_maps/names_indices', 'rb') as f:
    indices, residues, bead_sort, names = pickle.load(f)

results_file = open(f'{save_path}/saxs_end_to_end_summary.txt', 'w')

def write_and_print(text):
    print(text)
    results_file.write(text + '\n')

all_dist_arrays = []
all_labels      = []
all_targets     = []
summary_rows    = []

for restraint in tqdm.tqdm(saxs_dmax, desc='Calculating Global Dmax'):
    label  = restraint['label']
    prot   = restraint['prot']
    res_s  = restraint['res_start']
    res_e  = restraint['res_end']
    target = restraint['target_distance']
    thresh = restraint['threshold']

    idx = [indices[prot][i] for i in range(len(indices[prot])) if res_s <= residues[prot][i] <= res_e]
    idx = np.unique(idx)
    coords_subset = m[prot][:, idx, :]

    # --- NEW GLOBAL LOGIC ---
    # For every frame, find the MAXIMUM distance between any two beads in the subset
    dist_per_frame = []
    for f in range(coords_subset.shape[0]):
        frame_xyz = coords_subset[f, :, :]
        # pdist calculates all pairwise distances; .max() is the Dmax for this frame
        dmax = pdist(frame_xyz).max()
        dist_per_frame.append(dmax)
    
    dist_per_frame = np.array(dist_per_frame)
    # -------------------------

    # Stats
    mean_d, median_d = np.mean(dist_per_frame), np.median(dist_per_frame)
    sd_d = np.std(dist_per_frame)
    p5, p95 = np.percentile(dist_per_frame, [5, 95])
    
    lower, upper = target - thresh, target + thresh
    pct_satisfied = 100.0 * np.sum((dist_per_frame >= lower) & (dist_per_frame <= upper)) / len(dist_per_frame)

    write_and_print(f'# {label}\nMean: {mean_d:.2f} Å | Target: {target} ± {thresh} Å | Sat: {pct_satisfied:.1f}%\n')

    all_dist_arrays.append(dist_per_frame)
    all_labels.append(label)
    all_targets.append(target)
    summary_rows.append({'label': label, 'mean': mean_d, 'median': median_d, 'sd': sd_d, 'p5': p5, 'p95': p95, 'target': target, 'pct_satisfied': pct_satisfied})

# Outlier removal (top 1%) for cleaner plotting
all_dist_arrays_filtered = [np.array(x)[np.array(x) < np.percentile(np.array(x), 99)] for x in all_dist_arrays]
signed_ds_filtered = [all_dist_arrays_filtered[i] - all_targets[i] for i in range(len(all_dist_arrays_filtered))]

# --- PLOTTING: SIGNED DEVIATION ---
fig1, ax1 = plt.subplots(figsize=(10, 6))
v1 = ax1.violinplot(signed_ds_filtered, showmeans=False, showextrema=False, positions=np.arange(len(all_labels)))
ax1.axhline(0, color='black', lw=1, linestyle='--')
for j, part in enumerate(v1['bodies']):
    part.set_facecolor('#AED6F1')
    part.set_edgecolor('#1A5276')
    part.set_alpha(0.7)
    for perc in np.percentile(signed_ds_filtered[j], [5, 50, 95]):
        ax1.plot([j - 0.4, j + 0.4], [perc, perc], color='red', linestyle=':', lw=1.5)
ax1.set_xticks(np.arange(len(all_labels)))
ax1.set_xticklabels([l.replace('_', '\n') for l in all_labels])
ax1.set_ylabel('Model Dmax - SAXS Dmax (Å)')
ax1.set_title('Signed Deviation (Global Distance)')
plt.tight_layout()
plt.savefig(f'{save_path}/signed_saxs_end_to_end_validation.png')
plt.close()

# --- PLOTTING: ABSOLUTE DISTANCE ---
fig2, ax2 = plt.subplots(figsize=(10, 6))
v2 = ax2.violinplot(all_dist_arrays_filtered, showmeans=False, showextrema=False, positions=np.arange(len(all_labels)))
for j, part in enumerate(v2['bodies']):
    part.set_facecolor('#AED6F1')
    part.set_edgecolor('#1A5276')
    part.set_alpha(0.7)
    for perc in np.percentile(all_dist_arrays_filtered[j], [5, 50, 95]):
        ax2.plot([j - 0.4, j + 0.4], [perc, perc], color='red', linestyle=':', lw=1.5)
ax2.set_xticks(np.arange(len(all_labels)))
ax2.set_xticklabels([l.replace('_', '\n') for l in all_labels])
ax2.set_ylabel('Absolute Distance (Å)')
ax2.set_title('Absolute Dmax Distribution')
plt.tight_layout()
plt.savefig(f'{save_path}/saxs_end_to_end_validation.png')
plt.close()

results_file.close()
print(f'Validation complete. Files saved in {save_path}/')
