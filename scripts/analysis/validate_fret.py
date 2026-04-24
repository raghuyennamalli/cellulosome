import numpy as np
import matplotlib.pyplot as plt
import argparse
import pickle
import tqdm
import os

parser = argparse.ArgumentParser(
    description="FRET Distance Validation",
)
parser.add_argument("--rmf", type=str, required=True, help="Path to the combined RMF file.")
parser.add_argument("--output_dir", type=str, required=True, help="Path to the parent validation directory.")
args = parser.parse_args()

combined_rmf = args.rmf
save_path    = os.path.join(args.output_dir, "fret_validation")
os.makedirs(save_path, exist_ok=True)

distance_restraint_list_cellulosome = [
    {
        "prot1": "cipa", "residue1_start": 1475,  "residue1_end": 1475,
        "prot2": "cipa", "residue2_start": 1639, "residue2_end": 1639,
        "target_distance": 76.5, "threshold": 10.0,
        "label": "coh8 95 - coh9 260"
    },
    {
        "prot1": "cipa", "residue1_start": 1475, "residue1_end": 1475,
        "prot2": "cipa", "residue2_start": 1563, "residue2_end": 1563,
        "target_distance": 73.8, "threshold": 10.0,
        "label": "coh8 95 - coh9 183"
    }
]

# load pickle files
with open('./extracted_xyzr/saved_data', 'rb') as f:
    m, p = pickle.load(f)
with open('./contact_maps/names_indices', 'rb') as f:
    indices, residues, bead_sort, names = pickle.load(f)

print('Loaded pickle files!')
print('Molecules found:', list(m.keys()))
print('cipa shape:', m['cipa'].shape)

# open results file
results_file = open(f'{save_path}/distance_restraint_not_used_results.txt', 'w')

def write_and_print(text):
    print(text)
    results_file.write(text + '\n')

# store all distances for combined violin plot
all_dist_arrays = []
all_labels      = []
all_targets     = []
all_lowers      = []
all_uppers      = []
summary_rows    = []

for restraint in tqdm.tqdm(distance_restraint_list_cellulosome, desc='Distance Restraints'):
    label  = restraint['label']
    prot1  = restraint['prot1']
    prot2  = restraint['prot2']
    res1_s = restraint['residue1_start']
    res1_e = restraint['residue1_end']
    res2_s = restraint['residue2_start']
    res2_e = restraint['residue2_end']
    target = restraint['target_distance']
    thresh = restraint['threshold']
    lower  = target - thresh
    upper  = target + thresh

    # get bead indices exactly like original script
    local_idx1 = np.unique([
        indices[prot1][i]
        for i in range(len(indices[prot1]))
        if res1_s <= residues[prot1][i] <= res1_e
    ])
    local_idx2 = np.unique([
        indices[prot2][i]
        for i in range(len(indices[prot2]))
        if res2_s <= residues[prot2][i] <= res2_e
    ])

    assert len(local_idx1) > 0, f'No beads found for {prot1} res {res1_s}-{res1_e}'
    assert len(local_idx2) > 0, f'No beads found for {prot2} res {res2_s}-{res2_e}'

    # get coords for all frames at once - same as original script
    coords1 = m[prot1][:, local_idx1, :]
    coords2 = m[prot2][:, local_idx2, :]

    # center to center distance per frame
    center1 = coords1.mean(axis=1)
    center2 = coords2.mean(axis=1)
    dist_per_frame = np.linalg.norm(center1 - center2, axis=1)
    
    # statistics
    mean_d   = np.mean(dist_per_frame)
    median_d = np.median(dist_per_frame)
    sd_d     = np.std(dist_per_frame)
    p5, p95  = np.percentile(dist_per_frame, [5, 95])

    # satisfaction check
    satisfied     = np.sum((dist_per_frame >= lower) & (dist_per_frame <= upper))
    pct_satisfied = 100.0 * satisfied / len(dist_per_frame)

    write_and_print(f'# {label}')
    write_and_print(f'mean: {mean_d:.2f} Å, median: {median_d:.2f} Å, sd: {sd_d:.2f} Å')
    write_and_print(f'p5: {p5:.2f} Å, p95: {p95:.2f} Å')
    write_and_print(f'target: {target:.2f} ± {thresh:.2f} Å')
    write_and_print(f'% satisfied: {pct_satisfied:.1f}%')
    write_and_print('')

    all_dist_arrays.append(dist_per_frame)
    all_labels.append(label)
    all_targets.append(target)
    all_lowers.append(lower)
    all_uppers.append(upper)
    summary_rows.append({
        'label': label, 'mean': mean_d, 'median': median_d,
        'sd': sd_d, 'p5': p5, 'p95': p95,
        'target': target, 'lower': lower, 'upper': upper,
        'pct_satisfied': pct_satisfied,
    })

#outlier removal
all_dist_arrays_filtered = [np.array(x)[np.array(x) < np.percentile(np.array(x), 99)] for x in all_dist_arrays]
signed_ds_filtered = [all_dist_arrays_filtered[i] - all_targets[i] for i in range(len(all_dist_arrays_filtered))]

# 1. Create the Signed Distance Plot (Deviation from Target)
fig_signed, ax_s = plt.subplots(figsize=(10, 6))
vplot_s = ax_s.violinplot([np.array(x) for x in signed_ds_filtered],
                          showmeans=False, showextrema=False,
                          positions=np.arange(len(signed_ds_filtered)))

ax_s.axhline(0, color='black', lw=1, linestyle='--') 

for j, part in enumerate(vplot_s['bodies']):
    part.set_facecolor('#AED6F1')
    part.set_edgecolor('#1A5276')
    part.set_alpha(0.7)
    # Extended Red Percentiles (0.4 width)
    for perc in np.percentile(np.array(signed_ds_filtered[j]), [5, 50, 95]):
        ax_s.plot([j - 0.4, j + 0.4], [perc, perc], color='red', linestyle=':', lw=1.5)

ax_s.set_xticks(np.arange(len(all_labels)))
ax_s.set_xticklabels([l.replace('_', '\n') for l in all_labels], fontsize=9)
ax_s.set_ylabel('Measured - Target (Å)')
ax_s.set_title('FRET Deviation Summary (Signed)', fontsize=14)
plt.tight_layout()
plt.savefig(f'{save_path}/fret_signed_validation.png', dpi=150)
plt.close(fig_signed)

# 2. Create the Absolute Distance Plot
fig_abs, ax_a = plt.subplots(figsize=(10, 6))
vplot_a = ax_a.violinplot(all_dist_arrays_filtered,
                          showmeans=False, showextrema=False,
                          positions=np.arange(len(all_dist_arrays_filtered)))

for j, part in enumerate(vplot_a['bodies']):
    part.set_facecolor('#AED6F1')
    part.set_edgecolor('#1A5276')
    part.set_alpha(0.7)
    # Extended Red Percentiles (0.4 width)
    for perc in np.percentile(np.array(all_dist_arrays_filtered[j]), [5, 50, 95]):
        ax_a.plot([j - 0.4, j + 0.4], [perc, perc], color='red', linestyle=':', lw=1.5)

ax_a.set_xticks(np.arange(len(all_labels)))
ax_a.set_xticklabels([l.replace('_', '\n') for l in all_labels], fontsize=9)
ax_a.set_ylabel('Absolute Distance (Å)')
ax_a.set_title('FRET Absolute Distances', fontsize=14)
plt.tight_layout()
plt.savefig(f'{save_path}/fret_absolute_validation.png', dpi=150)
plt.close(fig_abs)

# summary table
write_and_print('# Summary')
header = (f"{'Label':<35} {'Mean':>8} {'Median':>8} {'SD':>8} "
          f"{'P5':>8} {'P95':>8} {'Target':>8} "
          f"{'Lower':>8} {'Upper':>8} {'%Sat':>8}")
write_and_print(header)
for r in summary_rows:
    write_and_print(
        f"{r['label']:<35} {r['mean']:>8.2f} {r['median']:>8.2f} "
        f"{r['sd']:>8.2f} {r['p5']:>8.2f} {r['p95']:>8.2f} "
        f"{r['target']:>8.2f} {r['lower']:>8.2f} {r['upper']:>8.2f} "
        f"{r['pct_satisfied']:>8.1f}")

results_file.close()
print('Distance Restraint Validation Complete!')
print(f'Results saved to {save_path}/')
