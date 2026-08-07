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
    description="Cela_SAXS_compare_pdist"
)
parser.add_argument("--h5", type=str, required=True, help="Path to the extracted XYZR HDF5 file.")
parser.add_argument("--output_dir", type=str, required=True, help="Path to the parent validation directory.")
args = parser.parse_args()

h5_file = args.h5
save_path    = os.path.join(args.output_dir, "Cela_SAXS_compare_pdist")
os.makedirs(save_path, exist_ok=True)
plt.rcParams['savefig.dpi'] = 600

saxs_dmax = [
    {"prot": "cela", "copy": 0, "res_start": 33, "res_end": 477, "target_distance": 118, "threshold": 3.0, "label": "cel8a_1"},
    {"prot": "cela", "copy": 1, "res_start": 33, "res_end": 477, "target_distance": 118, "threshold": 3.0, "label": "cel8a_2"},
    {"prot": "cela", "copy": 2, "res_start": 33, "res_end": 477, "target_distance": 118, "threshold": 3.0, "label": "cel8a_3"},
    {"prot": "cela", "copy": 3, "res_start": 33, "res_end": 477, "target_distance": 118, "threshold": 3.0, "label": "cel8a_4"},
    {"prot": "cela", "copy": 4, "res_start": 33, "res_end": 477, "target_distance": 118, "threshold": 3.0, "label": "cel8a_5"},
    {"prot": "cela", "copy": 5, "res_start": 33, "res_end": 477, "target_distance": 118, "threshold": 3.0, "label": "cel8a_6"},
    {"prot": "cela", "copy": 6, "res_start": 33, "res_end": 477, "target_distance": 118, "threshold": 3.0, "label": "cel8a_7"},
    {"prot": "cela", "copy": 7, "res_start": 33, "res_end": 477, "target_distance": 118, "threshold": 3.0, "label": "cel8a_8"},
    {"prot": "cela", "copy": 8, "res_start": 33, "res_end": 477, "target_distance": 118, "threshold": 3.0, "label": "cel8a_9"},
]

m_temp = defaultdict(list)
indices = defaultdict(list)
residues = defaultdict(list)
bead_copies = defaultdict(list) 

with h5py.File(h5_file, "r") as f:
    def parse_key_properties(key_string):
        parts = key_string.split("_")
        prot = parts[0]
        copy = int(parts[1]) if (len(parts) > 1 and parts[1].isdigit()) else 0
        res_part = parts[-1]
        match = re.search(r'([0-9]+)', res_part)
        res_start = int(match.group(1)) if match else 0
        return prot, copy, res_start

    sorted_keys = sorted(f.keys(), key=parse_key_properties)
    bead_counter_per_mol = defaultdict(int)

    for bead_key in sorted_keys:
        parts = bead_key.split("_")
        prot_name = parts[0]     
        res_range = parts[-1]    
        
        # Pull the (Frames, 3) raw coordinate matrix directly
        prot_name, copy_num, _ = parse_key_properties(bead_key)
        res_range = bead_key.split("_")[-1]
        
        xyz_coords = f[bead_key][:, 0:3] 
        m_temp[prot_name].append(xyz_coords)
            
        current_bead_idx = bead_counter_per_mol[prot_name]
        bead_counter_per_mol[prot_name] += 1
        bead_copies[prot_name].append(copy_num) 

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
for prot_name in bead_copies:
    bead_copies[prot_name] = np.array(bead_copies[prot_name])

results_file = open(f'{save_path}/Cela_SAXS_compare_pdist.txt', 'w')

def write_and_print(text):
    print(text)
    results_file.write(text + '\n')

all_dist_arrays = []
all_labels      = []
all_targets     = []
all_thresholds  = []
summary_rows    = []

for restraint in saxs_dmax:
    
    label = restraint['label']
    prot = restraint['prot']
    copy_id = restraint['copy'] 
    res_s = restraint['res_start']
    res_e = restraint['res_end']
    target = restraint['target_distance']
    thresh = restraint['threshold']
    
    idx = [indices[prot][i] for i in range(len(indices[prot])) if res_s <= residues[prot][i] <= res_e]
    idx = np.unique(idx)
        
    idx = np.array([i for i in idx if bead_copies[prot][i] == copy_id])
    if len(idx) == 0: 
        continue
            
    coords_subset = m[prot][:, idx, :]

    dist_per_frame = []
    for f_idx in range(coords_subset.shape[0]):
        frame_xyz = coords_subset[f_idx, :, :]
        dmax = pdist(frame_xyz).max()
        dist_per_frame.append(dmax)
    dist_per_frame = np.array(dist_per_frame)

    mean_d, median_d = np.mean(dist_per_frame), np.median(dist_per_frame)
    sd_d = np.std(dist_per_frame)
    p5, p95 = np.percentile(dist_per_frame, [5, 95])
        
    lower, upper = target - thresh, target + thresh
    pct_satisfied = 100.0 * np.sum((dist_per_frame >= lower) & (dist_per_frame <= upper)) / len(dist_per_frame)

    write_and_print(f'# {label}\nMean: {mean_d:.2f} Å | Median: {median_d:.2f} Å | Target: {target} ± {thresh} Å | Sat: {pct_satisfied:.1f}%\n')

    _lower_bound = target - thresh
    _upper_bound = target + thresh
    _total_f = len(dist_per_frame) if len(dist_per_frame) > 0 else 1
    
    _pct_below = 100.0 * np.sum(dist_per_frame < _lower_bound) / _total_f
    _pct_above = 100.0 * np.sum(dist_per_frame > _upper_bound) / _total_f
    
    write_and_print(
        f'-> [{label} Distribution Summary]\n'
        f'   Frames below target range (<{_lower_bound:.1f} Å): {_pct_below:.1f}%\n'
        f'   Frames above target range (>{_upper_bound:.1f} Å): {_pct_above:.1f}%\n'
    )
    
    all_dist_arrays.append(dist_per_frame)
    all_labels.append(label)
    all_targets.append(target)
    all_thresholds.append(thresh)
    summary_rows.append({'label': label, 'mean': mean_d, 'median': median_d, 'sd': sd_d, 'p5': p5, 'p95': p95, 'target': target, 'pct_satisfied': pct_satisfied})

all_dist_arrays_filtered = [np.array(x)[np.array(x) < np.percentile(np.array(x), 99)] for x in all_dist_arrays]
signed_ds_filtered = [all_dist_arrays_filtered[i] - all_targets[i] for i in range(len(all_dist_arrays_filtered))]

all_dist_filtered = [np.array(x) for x in all_dist_arrays]
if all_dist_filtered:
    fig, ax = plt.subplots(figsize=(10, 6))
    vplot = ax.violinplot(all_dist_filtered, showmeans=False, showextrema=False, positions=np.arange(len(all_labels)))
    
    for j, part in enumerate(vplot['bodies']):
        part.set_facecolor('#AED6F1'); part.set_edgecolor('#1A5276'); part.set_alpha(1.0)
        
        median_val = np.percentile(all_dist_filtered[j], 50)
        #ax.plot([j - 0.25, j + 0.25], [median_val, median_val], color='black', linestyle='-', lw=2.5, zorder=100)
        
        tgt = all_targets[j]
        thresh_val = all_thresholds[j] 
        ax.plot([j - 0.4, j + 0.4], [tgt - thresh_val, tgt - thresh_val], color='black', linestyle=':', lw=1.5, zorder=50)
        ax.plot([j - 0.4, j + 0.4], [tgt + thresh_val, tgt + thresh_val], color='black', linestyle=':', lw=1.5, zorder=50)
            
    raw_min = min(d.min() for d in all_dist_filtered)
    raw_max = max(d.max() for d in all_dist_filtered)
    padding = (raw_max - raw_min) * 0.15 if (raw_max - raw_min) > 0 else 15
    ymin = raw_min - padding
    ymax = raw_max + padding

    ax.set_xticks(np.arange(len(all_labels)))
    ax.set_xticklabels([l.replace('_', '\n') for l in all_labels], fontsize=9)
    ax.set_ylabel('Distance (Å)')
    ax.set_title('Dmax Distribution (Per Copy)')
    ax.set_ylim(ymin, ymax)
    plt.tight_layout()
    plt.savefig(f'{save_path}/Cela_SAXS_compare_pdist.png')
    plt.savefig(f'{save_path}/Cela_SAXS_compare_pdist.svg')
    plt.close()
results_file.close()
print(f'Validation complete. Files saved in {save_path}/')
