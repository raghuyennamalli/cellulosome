import numpy as np
import matplotlib.pyplot as plt
import h5py
import re
from collections import defaultdict
import tqdm
import os
import argparse

parser = argparse.ArgumentParser(description="Distance Restraint Validation")
parser.add_argument("--h5", type=str, required=True, help="Path to the extracted XYZR HDF5 file.")
parser.add_argument("--output_dir", type=str, required=True, help="Path to the parent validation directory.")
args = parser.parse_args()

h5_file = args.h5
save_path    = os.path.join(args.output_dir, "Distance_restraint_validation")
os.makedirs(save_path, exist_ok=True)
plt.rcParams['savefig.dpi'] = 600

# specific residue definitions
distance_restraint_list_cellulosome = [
    {"prot1": "cipa", "copy1": 0, "residue1_start": 967,  "residue1_end": 967,  "prot2": "cipa", "copy2": 0, "residue2_start": 1132, "residue2_end": 1132, "target_distance": 58.743, "threshold": 5.0, "label": "coh5 967 - coh6 1132"},
    {"prot1": "cipa", "copy1": 0, "residue1_start": 1132, "residue1_end": 1132, "prot2": "cipa", "copy2": 0, "residue2_start": 1297, "residue2_end": 1297, "target_distance": 58.743, "threshold": 5.0, "label": "coh6 1132 - coh7 1297"},
    {"prot1": "cipa", "copy1": 0, "residue1_start": 1297, "residue1_end": 1297, "prot2": "cipa", "copy2": 0, "residue2_start": 1462, "residue2_end": 1462, "target_distance": 58.743, "threshold": 5.0, "label": "coh7 1297 - coh8 1462"},
    {"prot1": "cipa", "copy1": 0, "residue1_start": 1462, "residue1_end": 1462, "prot2": "cipa", "copy2": 0, "residue2_start": 1626, "residue2_end": 1626, "target_distance": 58.743, "threshold": 5.0, "label": "coh8 1462 - coh9 1626"}
]

 # Create temporary lists to safely collect arrays before turning them into a NumPy matrix

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
results_file = open(f'{save_path}/Distance_restraint_results.txt', 'w')
def write_and_print(text):
    print(text)
    results_file.write(text + '\n')

all_dist_arrays, all_labels, summary_stats = [], [], []

for restraint in tqdm.tqdm(distance_restraint_list_cellulosome, desc='Processing'):
    label, prot1, prot2 = restraint['label'], restraint['prot1'], restraint['prot2']
    c1, c2 = restraint['copy1'], restraint['copy2']
    r1s, r1e, r2s, r2e = restraint['residue1_start'], restraint['residue1_end'], restraint['residue2_start'], restraint['residue2_end']
    target, thresh = restraint['target_distance'], restraint['threshold']
    
    idx1 = np.unique([indices[prot1][i] for i in range(len(indices[prot1])) if r1s <= residues[prot1][i] <= r1e])
    idx2 = np.unique([indices[prot2][i] for i in range(len(indices[prot2])) if r2s <= residues[prot2][i] <= r2e])
    
    idx1 = np.array([idx for idx in idx1 if bead_copies[prot1][idx] == c1])
    idx2 = np.array([idx for idx in idx2 if bead_copies[prot2][idx] == c2])
    
    dist = np.linalg.norm(m[prot1][:, idx1, :].mean(axis=1) - m[prot2][:, idx2, :].mean(axis=1), axis=1)
    
    lower, upper = target - thresh, target + thresh
    satisfied = np.sum((dist >= lower) & (dist <= upper))
    pct_satisfied = 100.0 * satisfied / len(dist)
    
    all_dist_arrays.append(dist)
    all_labels.append(label)
    summary_stats.append({
        'label': label, 'mean': np.mean(dist), 'median': np.median(dist), 
        'target': target, 'pct_satisfied': pct_satisfied
    })

all_dist_filtered = [np.array(x)[np.array(x) < np.percentile(np.array(x), 99)] for x in all_dist_arrays]
all_targets = [s['target'] for s in summary_stats]
signed_ds_filtered = [all_dist_filtered[i] - all_targets[i] for i in range(len(all_dist_filtered))]

write_and_print(f"{'Label':<35} {'Mean':>10} {'Median':>10} {'Target':>10} {'%Satisfied':>12}")
write_and_print("-" * 81)
for s in summary_stats:
    write_and_print(f"{s['label']:<35} {s['mean']:>10.2f} {s['median']:>10.2f} {s['target']:>10.2f} {s['pct_satisfied']:>11.1f}%")

def create_violin_plot(data_list, label_list, filename, title, ymin, ymax, targets_list, thresholds_list):
    fig, ax = plt.subplots(figsize=(8, 6))
    vplot = ax.violinplot([np.array(x) for x in data_list], showmeans=False, showextrema=False, positions=np.arange(len(data_list)))
    
    for j, part in enumerate(vplot['bodies']):
        part.set_facecolor('#AED6F1'); part.set_edgecolor('#1A5276'); part.set_alpha(1.0)
        for perc in np.percentile(np.array(data_list[j]), [5, 50, 95]):
            ax.plot([j - 0.25, j + 0.25], [perc, perc], color='red', linestyle=':', lw=2, zorder=100)
        t_val = targets_list[j]
        #thresh_val = thresholds_list[j]
        #ax.plot([j - 0.4, j + 0.4], [t_val - thresh_val, t_val - thresh_val], color='black', linestyle=':', lw=1.5, zorder=50)
        #ax.plot([j - 0.4, j + 0.4], [t_val + thresh_val, t_val + thresh_val], color='black', linestyle=':', lw=1.5, zorder=50)
            
    ax.set_xticks(np.arange(len(label_list)))
    ax.set_xticklabels([l.replace('_', '\n') for l in label_list], fontsize=9)
    ax.set_ylabel('Center-to-center distance (Å)')
    ax.set_title(title)
    ax.set_ylim(ymin, ymax)
    plt.tight_layout()
    plt.savefig(f'{save_path}/{filename}')
    plt.close()

all_thresholds = [s['threshold'] if 'threshold' in s else 5.0 for s in distance_restraint_list_cellulosome]

coh_data = [d for d, l in zip(all_dist_filtered, all_labels) if "coh" in l]
coh_lbls = [l for l in all_labels if "coh" in l]
coh_tgts = [t for t, l in zip(all_targets, all_labels) if "coh" in l]
coh_thrs = [t for t, l in zip(all_thresholds, all_labels) if "coh" in l]

if coh_data:
    create_violin_plot(coh_data, coh_lbls, 'coh_absolute_restraints.png', 'Cohesin Distance Restraints', 50, 70, coh_tgts, coh_thrs)
    create_violin_plot(coh_data, coh_lbls, 'coh_absolute_restraints.svg', 'Cohesin Distance Restraints', 50, 70, coh_tgts, coh_thrs)

fig_s, ax_s = plt.subplots(figsize=(8, 6))
vplot_s = ax_s.violinplot([np.array(x) for x in signed_ds_filtered],
                          showmeans=False, showextrema=False,
                          positions=np.arange(len(signed_ds_filtered)))
ax_s.axhline(0, color='black', lw=1, linestyle='--')
ax_s.axhline(-5.0, color='black', lw=1.5, linestyle=':')
ax_s.axhline(5.0, color='black', lw=1.5, linestyle=':')
for j, part in enumerate(vplot_s['bodies']):
    part.set_facecolor('#AED6F1')
    part.set_edgecolor('#1A5276')
    part.set_alpha(1.0)
    for perc in np.percentile(np.array(signed_ds_filtered[j]), [5, 50, 95]):
        ax_s.plot([j - 0.25, j + 0.25], [perc, perc], color='red', linestyle=':', lw=2, zorder=100)
ax_s.set_xticks(np.arange(len(all_labels)))
ax_s.set_xticklabels([l.replace('_', '\n').replace(' - ', '\n') for l in all_labels], fontsize=9)
ax_s.set_ylabel('Measured - Target (Å)', fontweight='normal')
ax_s.set_title('EM Derived Distance Restraints (Signed Deviation)', fontsize=12)
plt.tight_layout()
plt.savefig(f'{save_path}/em_derived_restraints_signed.svg')
plt.savefig(f'{save_path}/em_derived_restraints_signed.png')
plt.close(fig_s)
results_file.close()
print(f"Success! Images and '{save_path}/distance_restraint_results.txt' generated.")

