
# Integrative modeling of _Acetivibrio thermocellus_ cellulosome

This repository consists of input files and scripts used for the integrative modelling of the _Acetivibrio thermocellus_ cellulosome. The major cluster center model and localization probability density maps were also added. The cellulosome model is based on the AlphaFold predicted structures, crystal structures, cryo-EM, SAXS, connectivity and stereochemistry. The modeling was performed using [IMP](https://integrativemodeling.org) (*Integrative Modeling Platform*).

## Directory structure
1. [inputs](inputs/) : Contains the subdirectories for the input data used for the modeling.
2. [scripts](scripts/) : Contains all the scripts used for modeling, analysis and validation of the models.
3. [results](results/) : Contains the major cluster center model and the localization probability densities of the top cluster of the subcomplexes.

## Protocol
### Preprocessing
[af_pipeline](https://github.com/isblab/af_pipeline): Used to rank and extract rigid bodies from the AlphaFold3 predicted structures

### Sampling
To perform sampling:
```
./run_modeling.sh <output_dir> <number_of_runs> <number_of_steps>
```

where, \
`output_dir` Path to the output directory \
`number_of_runs` Total number of sampling runs to perform\
`number_of_steps` Total number of steps

### Analysis

Analysis was performed with the scripts in [IMP_Toolbox](https://github.com/isblab/IMP_Toolbox/tree/main/IMP_Toolbox/analysis).

`run_analysis_trajectories.py`: Select good scoring models \
`variable_filter.py`: Filter the good scoring models as per the model cap (30,000 in this case) \ 
`run_extract_models.py`: Extract the filtered good scoring models \
`exhaust.py`: Analyse the sampling convergence \
`extract_sampcon.py`: Extract the models in the major cluster center model \
`rmf_to_xyzr.py`: Convert and save RMF file to XYZR format\
`correlate_clsiter_sample_densities.py`:Correlate Sample A and Sample B LPD maps\
`interaction_map.py`: Parse XYZR data and calculate the contacts between the domains\
`align_pdb_to_ccm.py`: Fit the atomic structures used for modelling onto the major cluster center model\
[PrISM](https://doi.org/10.1093/bioinformatics/btac400): Obtain Domain-wise precision of the major cluster
`end_to_end_analysis.py` script was used to automate the analysis

`density_domains_linkers.txt` was used to generate LPDs.

### Validation
To run the cryo-EM validation: 
```
python fit_to_em_data.py --output_dir path/to/output/folder
``` 

To run the validation scripts (`fit_to_saxs_end_to_end.py`, `fit_to_distance_restraint.py`, `validate_fret.py`):
```
python <file_name>.py --h5 path/to/h5/file --output_dir path/to/output/folder
```

### Results

For each of the simulations, the following files are in the [results](results/) directory
* `cluster_center_model.rmf3` : representative bead model of the major cluster
* `LPD_*.mrc` : Localization probability density maps
* `aligned_*.pdb` : Atomic structures fitted to the major cluster center model

## Addtional information

**Author(s):** Nisha Nandhini Shankar, Omkar Golatkar, Shruthi Viswanath, U. Venkatasubramanian, Ragothaman M. Yennamalli\
**Date:**  21st August 2026\
**License**: [GPLv3](LICENSE.txt)\
**Publications:** Integrative structural modeling reveals the global architecture and conformational heterogeneity of _Acetivibrio thermocellus_ cellulosome
