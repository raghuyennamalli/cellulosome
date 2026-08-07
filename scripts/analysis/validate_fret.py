import numpy as np
import matplotlib.pyplot as plt
import argparse
import h5py
import re
from collections import defaultdict
import tqdm
import os

parser = argparse.ArgumentParser(
    description="FRET Distance Validation",
)
parser.add_argument("--h5", type=str, required=True, help="Path to the extracted XYZR HDF5 file.")
parser.add_argument("--output_dir", type=str, required=True, help="Path to the parent validation directory.")
args = parser.parse_args()

h5_file = args.h5
save_path    = os.path.join(args.output_dir, "fret_validation")
os.makedirs(save_path, exist_ok=True)

plt.rcParams['savefig.dpi'] = 600


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
results_file = open(f'{save_path}/distance_restraint_not_used_results.txt', 'w')

def write_and_print(text):
    print(text)
    results_file.write(text + '\n')

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

    coords1 = m[prot1][:, local_idx1, :]
    coords2 = m[prot2][:, local_idx2, :]

    center1 = coords1.mean(axis=1)
    center2 = coords2.mean(axis=1)
    dist_per_frame = np.linalg.norm(center1 - center2, axis=1)
    
    mean_d   = np.mean(dist_per_frame)
    median_d = np.median(dist_per_frame)
    sd_d     = np.std(dist_per_frame)
    p5, p95  = np.percentile(dist_per_frame, [5, 95])

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

all_dist_arrays_filtered = [np.array(x)[np.array(x) < np.percentile(np.array(x), 99)] for x in all_dist_arrays]
signed_ds_filtered = [all_dist_arrays_filtered[i] - all_targets[i] for i in range(len(all_dist_arrays_filtered))]

fig_signed, ax_s = plt.subplots(figsize=(10, 6))
vplot_s = ax_s.violinplot([np.array(x) for x in signed_ds_filtered],
                          showmeans=False, showextrema=False,
                          positions=np.arange(len(signed_ds_filtered)))

ax_s.axhline(0, color='black', lw=1, linestyle='--') 

for j, part in enumerate(vplot_s['bodies']):
    part.set_facecolor('#AED6F1')
    part.set_edgecolor('#1A5276')
    part.set_alpha(1.0)
    for perc in np.percentile(np.array(signed_ds_filtered[j]), [5, 50, 95]):
        ax_s.plot([j - 0.4, j + 0.4], [perc, perc], color='red', linestyle=':', lw=1.5)
    
    thresh = distance_restraint_list_cellulosome[j]['threshold']
    ax_s.axhline(thresh, color='black', linestyle='--', alpha=0.7, label='Satisfaction Bounds' if j==0 else "")
    ax_s.axhline(-thresh, color='black', linestyle='--', alpha=0.7)


ax_s.set_xticks(np.arange(len(all_labels)))
ax_s.set_xticklabels([l.replace('_', '\n') for l in all_labels], fontsize=9)
ax_s.set_ylabel('Measured - Target (Å)')
ax_s.set_title('FRET Deviation Summary (Signed)', fontsize=14)
plt.tight_layout()
plt.savefig(f'{save_path}/fret_signed_validation.png', dpi=150)
plt.close(fig_signed)

fig_abs, ax_a = plt.subplots(figsize=(10, 6))
vplot_a = ax_a.violinplot(all_dist_arrays_filtered,
                          showmeans=False, showextrema=False,
                          positions=np.arange(len(all_dist_arrays_filtered)))

for j, part in enumerate(vplot_a['bodies']):
    part.set_facecolor('#AED6F1')
    part.set_edgecolor('#1A5276')
    part.set_alpha(1.0)
    for perc in np.percentile(np.array(all_dist_arrays_filtered[j]), [5, 50, 95]):
        ax_a.plot([j - 0.4, j + 0.4], [perc, perc], color='red', linestyle=':', lw=1.5)
    #low, upp = all_lowers[j], all_uppers[j]
    #ax_a.plot([j - 0.4, j + 0.4], [low, low], color='black', linestyle='--', alpha=0.7, label='Satisfaction Bounds' if j==0 else "")
    #ax_a.plot([j - 0.4, j + 0.4], [upp, upp], color='black', linestyle='--', alpha=0.7)

ax_a.set_xticks(np.arange(len(all_labels)))
ax_a.set_xticklabels([l.replace('_', '\n') for l in all_labels], fontsize=9)
ax_a.set_ylabel('Distance (Å)')
ax_a.set_title('FRET Absolute Distances', fontsize=14)
plt.tight_layout()
plt.savefig(f'{save_path}/fret_absolute_validation.png', dpi=150)
plt.close(fig_abs)

write_and_print('# Summary')
header = (f"{'Label':<35} {'Mean':>8} {'Median':>8} {'SD':>8} "
          f"{'P5':>8} {'P95':>8} {'Target':>8} "
          f"{'Lower':>8} {'Upper':>8} {'%Sat':>8}")
write_and_print(header)
write_and_print("-" * len(header))
for r in summary_rows:
    write_and_print(
        f"{r['label']:<35} {r['mean']:>8.2f} {r['median']:>8.2f} "
        f"{r['sd']:>8.2f} {r['p5']:>8.2f} {r['p95']:>8.2f} "
        f"{r['target']:>8.2f} {r['lower']:>8.2f} {r['upper']:>8.2f} "
        f"{r['pct_satisfied']:>8.1f}")

results_file.close()
print('Distance Restraint Validation Complete!')
print(f'Results saved to {save_path}/')
