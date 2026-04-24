import numpy as np
import matplotlib.pyplot as plt
from rmf_xyzr import get_all_data_xyzr
import pickle
import tqdm
import os
import argparse

parser = argparse.ArgumentParser(description="Distance Restraint Validation")
parser.add_argument(
    "--rmf",
    type=str,
    required=True,
    help="Path to the combined RMF file."
)
parser.add_argument(
    "--output_dir",
    type=str,
    required=True,
    help="Path to the parent validation directory."
)
args = parser.parse_args()

combined_rmf = args.rmf
save_path    = os.path.join(args.output_dir, "distance_restraint_validation")
os.makedirs(save_path, exist_ok=True)

# Define Restraints
distance_restraint_list_cellulosome = [
    {"prot1": "cipa", "residue1_start": 967,  "residue1_end": 967,  "prot2": "cipa", "residue2_start": 1132, "residue2_end": 1132, "target_distance": 58.743, "threshold": 5.0, "label": "coh5 967 - coh6 1132"},
    {"prot1": "cipa", "residue1_start": 1132, "residue1_end": 1132, "prot2": "cipa", "residue2_start": 1297, "residue2_end": 1297, "target_distance": 58.743, "threshold": 5.0, "label": "coh6 1132 - coh7 1297"},
    {"prot1": "cipa", "residue1_start": 1297, "residue1_end": 1297, "prot2": "cipa", "residue2_start": 1462, "residue2_end": 1462, "target_distance": 58.743, "threshold": 5.0, "label": "coh7 1297 - coh8 1462"},
    {"prot1": "cipa", "residue1_start": 1462, "residue1_end": 1462, "prot2": "cipa", "residue2_start": 1626, "residue2_end": 1626, "target_distance": 58.743, "threshold": 5.0, "label": "coh8 1462 - coh9 1626"},
]

#data extraction
get_all_data_xyzr([combined_rmf], 30)
with open('./extracted_xyzr/saved_data', 'rb') as f:
    m, p = pickle.load(f)
with open('./contact_maps/names_indices', 'rb') as f:
    indices, residues, bead_sort, names = pickle.load(f)

# Open results file for satisfaction percentages
results_file = open(f'{save_path}/distance_restraint_results.txt', 'w')
def write_and_print(text):
    print(text)
    results_file.write(text + '\n')

all_dist_arrays, all_labels, summary_stats = [], [], []

for restraint in tqdm.tqdm(distance_restraint_list_cellulosome, desc='Processing'):
    label, prot1, prot2 = restraint['label'], restraint['prot1'], restraint['prot2']
    r1s, r1e, r2s, r2e = restraint['residue1_start'], restraint['residue1_end'], restraint['residue2_start'], restraint['residue2_end']
    target, thresh = restraint['target_distance'], restraint['threshold']
    
    idx1 = np.unique([indices[prot1][i] for i in range(len(indices[prot1])) if r1s <= residues[prot1][i] <= r1e])
    idx2 = np.unique([indices[prot2][i] for i in range(len(indices[prot2])) if r2s <= residues[prot2][i] <= r2e])
    
    dist = np.linalg.norm(m[prot1][:, idx1, :].mean(axis=1) - m[prot2][:, idx2, :].mean(axis=1), axis=1)
    
    # Calculate Satisfaction Percentage
    lower, upper = target - thresh, target + thresh
    satisfied = np.sum((dist >= lower) & (dist <= upper))
    pct_satisfied = 100.0 * satisfied / len(dist)
    
    all_dist_arrays.append(dist)
    all_labels.append(label)
    summary_stats.append({
        'label': label, 'mean': np.mean(dist), 'median': np.median(dist), 
        'target': target, 'pct_satisfied': pct_satisfied
    })

# Remove extreme outliers (top 1% of distances) for cleaner plotting
all_dist_filtered = [np.array(x)[np.array(x) < np.percentile(np.array(x), 99)] for x in all_dist_arrays]
all_targets = [s['target'] for s in summary_stats]
signed_ds_filtered = [all_dist_filtered[i] - all_targets[i] for i in range(len(all_dist_filtered))]

# Write statistics to text file
write_and_print(f"{'Label':<35} {'Mean':>8} {'Median':>8} {'Target':>8} {'%Satisfied':>12}")
for s in summary_stats:
    write_and_print(f"{s['label']:<35} {s['mean']:>8.2f} {s['median']:>8.2f} {s['target']:>8.2f} {s['pct_satisfied']:>11.1f}%")

# Plotting Function
def create_violin_plot(data_list, label_list, filename, title, ymin, ymax):
    fig, ax = plt.subplots(figsize=(8, 6))
    vplot = ax.violinplot([np.array(x) for x in data_list], showmeans=False, showextrema=False, positions=np.arange(len(data_list)))
    
    for j, part in enumerate(vplot['bodies']):
        part.set_facecolor('#AED6F1'); part.set_edgecolor('#1A5276'); part.set_alpha(0.7)
        # Red Dotted Percentiles (5, 50, 95)
        for perc in np.percentile(np.array(data_list[j]), [5, 50, 95]):
            ax.plot([j - 0.25, j + 0.25], [perc, perc], color='red', linestyle=':', lw=2, zorder=100)
            
    ax.set_xticks(np.arange(len(label_list)))
    ax.set_xticklabels([l.replace('_', '\n') for l in label_list], fontsize=9)
    ax.set_ylabel('Center-to-center distance (Å)')
    ax.set_title(title)
    ax.set_ylim(ymin, ymax) # Custom Y-Limits applied here
    plt.tight_layout()
    plt.savefig(f'{save_path}/{filename}')
    plt.close()

# Absolute distance plots
COH_YMIN, COH_YMAX = 50, 70
create_violin_plot(all_dist_filtered, all_labels, 'em_derived_restraints.svg', 'EM Derived Distance Restraints', COH_YMIN, COH_YMAX)
create_violin_plot(all_dist_filtered, all_labels, 'em_derived_restraints.png', 'EM Derived Distance Restraints', COH_YMIN, COH_YMAX)

# Signed deviation plot
fig_s, ax_s = plt.subplots(figsize=(8, 6))
vplot_s = ax_s.violinplot([np.array(x) for x in signed_ds_filtered],
                          showmeans=False, showextrema=False,
                          positions=np.arange(len(signed_ds_filtered)))
ax_s.axhline(0, color='black', lw=1, linestyle='--')
for j, part in enumerate(vplot_s['bodies']):
    part.set_facecolor('#AED6F1')
    part.set_edgecolor('#1A5276')
    part.set_alpha(0.7)
    for perc in np.percentile(np.array(signed_ds_filtered[j]), [5, 50, 95]):
        ax_s.plot([j - 0.25, j + 0.25], [perc, perc], color='red', linestyle=':', lw=2, zorder=100)
ax_s.set_xticks(np.arange(len(all_labels)))
ax_s.set_xticklabels([l.replace('_', '\n') for l in all_labels], fontsize=9)
ax_s.set_ylabel('Measured - Target (Å)', fontweight='bold')
ax_s.set_title('EM Derived Distance Restraints (Signed Deviation)', fontsize=12)
plt.tight_layout()
plt.savefig(f'{save_path}/em_derived_restraints_signed.svg')
plt.savefig(f'{save_path}/em_derived_restraints_signed.png')
plt.close(fig_s)
results_file.close()
print(f"Success! Images and '{save_path}/distance_restraint_results.txt' generated.")

