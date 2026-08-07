import numpy as np
import matplotlib.pyplot as plt
import h5py
import re
from collections import defaultdict
import tqdm
import os
from scipy.spatial.distance import pdist
import argparse

parser = argparse.ArgumentParser(
    description="SAXS cohesin complexes comparison"
)
parser.add_argument("--h5", type=str, required=True, help="Path to the extracted XYZR HDF5 file.")
parser.add_argument("--output_dir", type=str, required=True, help="Path to the parent validation directory.")
args = parser.parse_args()

h5_file = args.h5
save_path    = os.path.join(args.output_dir, "saxs_cohesin_comparison")
os.makedirs(save_path, exist_ok=True)
plt.rcParams['savefig.dpi'] = 600

saxs_dmax = [
    {"prot": "cipa", "res_start": 29,    "res_end": 321, "target_distance": 207.0,  "threshold": 10.0, "label": "cohI-cohII"},
    {"prot": "cipa", "res_start": 184,    "res_end": 702, "target_distance": 232.0,  "threshold": 10.0, "label": "cohII-CBM-cohIII"}
]

m_temp = defaultdict(list)
indices = defaultdict(list)
residues = defaultdict(list)

with h5py.File(h5_file, "r") as f:
    def extract_start_res(key_string):
        res_part = key_string.split("_")[-1]
        match = re.search(r'([0-9]+)', res_part)
        return int(match.group(1)) if match else 0

    sorted_keys = sorted(f.keys(), key=lambda x: (x.split("_")[0], extract_start_res(x)))
    bead_counter_per_mol = defaultdict(int)

    for bead_key in sorted_keys:
        parts = bead_key.split("_")
        prot_name = parts[0]     
        res_range = parts[-1]    
        
        xyz_coords = f[bead_key][:, 0:3] 
        m_temp[prot_name].append(xyz_coords)
            
        current_bead_idx = bead_counter_per_mol[prot_name]
        bead_counter_per_mol[prot_name] += 1

        range_match = re.search(r'([0-9]+)-([0-9]+)', res_range)
        if range_match:
            n1, n2 = int(range_match.group(1)), int(range_match.group(2))
            res_list = list(range(n1, n2 + 1))
            residues[prot_name].extend(res_list)
            indices[prot_name].extend([current_bead_idx] * len(res_list))
        elif res_range.isdigit():  
            res_num = int(res_range)
            residues[prot_name].append(res_num)
            indices[prot_name].append(current_bead_idx)

m = {}
for prot_name, list_of_arrays in m_temp.items():
    m[prot_name] = np.stack(list_of_arrays, axis=1)

results_file = open(f'{save_path}/saxs_validation_summary.txt', 'w')

def write_and_print(text):
    print(text)
    results_file.write(text + '\n')

all_dist_arrays = []
all_labels      = []
all_targets     = []
summary_rows    = []

for restraint in tqdm.tqdm(saxs_dmax, desc='Calculating local Dmax'):
    label  = restraint['label']
    prot   = restraint['prot']
    res_s  = restraint['res_start']
    res_e  = restraint['res_end']
    target = restraint['target_distance']
    thresh = restraint['threshold']

    idx = [indices[prot][i] for i in range(len(indices[prot])) if res_s <= residues[prot][i] <= res_e]
    idx = np.unique(idx)
    coords_subset = m[prot][:, idx, :]

    dist_per_frame = []
    for f in range(coords_subset.shape[0]):
        frame_xyz = coords_subset[f, :, :]
        dmax = pdist(frame_xyz).max()
        dist_per_frame.append(dmax)
    
    dist_per_frame = np.array(dist_per_frame)

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

all_dist_arrays_filtered = [np.array(x)[np.array(x) < np.percentile(np.array(x), 99)] for x in all_dist_arrays]
signed_ds_filtered = [all_dist_arrays_filtered[i] - all_targets[i] for i in range(len(all_dist_arrays_filtered))]

fig1, ax1 = plt.subplots(figsize=(10, 6))
v1 = ax1.violinplot(signed_ds_filtered, showmeans=False, showextrema=False, positions=np.arange(len(all_labels)))
ax1.axhline(0, color='black', lw=1, linestyle='--')
for j, part in enumerate(v1['bodies']):
    part.set_facecolor('#AED6F1')
    part.set_edgecolor('#1A5276')
    part.set_alpha(1.0)
    for perc in np.percentile(signed_ds_filtered[j], [5, 50, 95]):
        ax1.plot([j - 0.4, j + 0.4], [perc, perc], color='red', linestyle=':', lw=1.5)
    #thresh = saxs_dmax[j]['threshold']
    #ax1.axhline(thresh, color='black', linestyle='--', alpha=0.6, label='Satisfaction Bounds' if j==0 else "")
    #ax1.axhline(-thresh, color='black', linestyle='--', alpha=0.6)
    
ax1.set_xticks(np.arange(len(all_labels)))
ax1.set_xticklabels([l.replace('_', '\n') for l in all_labels])
ax1.set_ylabel('Model Dmax - SAXS Dmax (Å)')
ax1.set_title('Signed Deviation')
plt.tight_layout()
plt.savefig(f'{save_path}/signed_saxs_validation.png')
plt.close()

fig2, ax2 = plt.subplots(figsize=(10, 6))
v2 = ax2.violinplot(all_dist_arrays_filtered, showmeans=False, showextrema=False, positions=np.arange(len(all_labels)))
for j, part in enumerate(v2['bodies']):
    part.set_facecolor('#AED6F1')
    part.set_edgecolor('#1A5276')
    part.set_alpha(0.7)
    for perc in np.percentile(all_dist_arrays_filtered[j], [5, 50, 95]):
        ax2.plot([j - 0.4, j + 0.4], [perc, perc], color='red', linestyle=':', lw=1.5)
    tgt = all_targets[j]
    #thresh = saxs_dmax[j]['threshold']
    #ax2.plot([j - 0.4, j + 0.4], [tgt - thresh, tgt - thresh], color='black', linestyle='--', alpha=0.7, label='Satisfaction Bounds' if j==0 else "")
    #ax2.plot([j - 0.4, j + 0.4], [tgt + thresh, tgt + thresh], color='black', linestyle='--', alpha=0.7)
    
ax2.set_xticks(np.arange(len(all_labels)))
ax2.set_xticklabels([l.replace('_', '\n') for l in all_labels])
ax2.set_ylabel('Distance (Å)')
ax2.set_title('Dmax Distribution')
plt.tight_layout()
plt.savefig(f'{save_path}/saxs_validation.png')
plt.close()

results_file.close()
print(f'Validation complete. Files saved in {save_path}/')
