import os
import numpy as np
import mrcfile
import argparse
from mrc_parser import calculate_with_external_grid_with_addition

# Paths
SEG_DIR = '/path/to/segmented_maps/segmented_normalized'
LPD_DIR = '/path/to/cluster.0'
parser = argparse.ArgumentParser(description="Whole Complex EM Validation")
parser.add_argument("--output_dir", type=str, required=True, help="Path to save validation results.")
args = parser.parse_args()

save_path = os.path.join(args.output_dir, "em_validation")
os.makedirs(save_path, exist_ok=True)

# 1. Load the full EMDB reference map
mrc_full = mrcfile.open(f'{SEG_DIR}/EMD_1823_coh_3_4_5_norm_inverted.mrc', 'r', permissive=True)

# 2. Define the new LPD filenames based on your list
lpd_filenames = [
    'LPD_Cohesin_3.mrc', 'LPD_Coh_3_coh_4_linker.mrc', 
    'LPD_Cohesin_4.mrc', 'LPD_Coh_4_coh_5_linker.mrc', 
    'LPD_Cohesin_5.mrc', 
    'LPD_Enzyme_3.mrc', 'LPD_Enzyme_3_linker.mrc', 'LPD_Dockerin_3.mrc', 
    'LPD_Enzyme_4.mrc', 'LPD_Enzyme_4_linker.mrc', 'LPD_Dockerin_4.mrc', 
    'LPD_Enzyme_5.mrc', 'LPD_Enzyme_5_linker.mrc', 'LPD_Dockerin_5.mrc'
    ]

# Helper to extract grid and values
def get_grid_and_values(mrc):
    xvals = mrc.voxel_size.x * np.arange(mrc.data.shape[2]) + mrc.header.origin.x
    yvals = mrc.voxel_size.y * np.arange(mrc.data.shape[1]) + mrc.header.origin.y
    zvals = mrc.voxel_size.z * np.arange(mrc.data.shape[0]) + mrc.header.origin.z
    values = mrc.data.transpose(2, 1, 0).copy()
    return (xvals, yvals, zvals), values

# Helper to save output MLPD_cela_3.mrcRC
def save_mrc(filepath, data, common_grid, voxel_spacing=5.0):
    with mrcfile.new(filepath, overwrite=True) as mrc_out:
        mrc_out.set_data(data.transpose(2, 1, 0).astype(np.float32))
        mrc_out.voxel_size = voxel_spacing
        mrc_out.header.origin.x = float(common_grid[0].min())
        mrc_out.header.origin.y = float(common_grid[1].min())
        mrc_out.header.origin.z = float(common_grid[2].min())

# Initialize lists for the comparison logic
list_of_g = []
list_of_v = []

# Load Full Map Grid/Values
g_f, v_f = get_grid_and_values(mrc_full)
list_of_g.append(g_f)
list_of_v.append(v_f)

# Load all LPD maps
for fname in lpd_filenames:
    path = os.path.join(LPD_DIR, fname)
    with mrcfile.open(path, 'r', permissive=True) as mrc:
        g, v = get_grid_and_values(mrc)
        list_of_g.append(g)
        list_of_v.append(v)

print(f'Loaded reference map and {len(lpd_filenames)} LPD components.')

# Whole Complex Comparison
results_path = os.path.join(save_path, 'em_validation_results_whole_complex.txt')
with open(results_path, 'w') as f:
    def write_and_print(text):
        print(text)
        f.write(text + '\n')

    write_and_print('# Whole Complex Comparison (Cela 0-8 + Cipa 0)')
    
    # This function interpolates all maps onto a common grid and sums the LPDs 
    # to compare against the first map in the list (the full reference).
    all_pts, nonzero_pts, summed_map, common_grid = calculate_with_external_grid_with_addition(
        list_of_g, list_of_v, 'Whole Complex', voxel_spacing=5.0)

    write_and_print(f'all points -> Overlap: {all_pts[0]:.4f}, C0: {all_pts[1]:.4f}, Cam: {all_pts[2]:.4f}, Pearson: {all_pts[3][0]:.4f}, Spearman: {all_pts[3][1]:.4f}')
    write_and_print(f'non-zero points -> Overlap: {nonzero_pts[0]:.4f}, C0: {nonzero_pts[1]:.4f}, Cam: {nonzero_pts[2]:.4f}, Pearson: {nonzero_pts[3][0]:.4f}, Spearman: {nonzero_pts[3][1]:.4f}')
    
    save_mrc(f'{save_path}/whole_complex_lpd_combined.mrc', summed_map, common_grid)

print(f'Done! Results saved to {save_path}')
