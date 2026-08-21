
# Integrative modeling of _Acetivibrio thermocellus_ cellulosome

This repository consists of input files and scripts used for the integrative modelling of the _Acetivibrio thermocellus_ cellulosome. The dominant cluster center model and localization probability density maps were also added. The cellulosome model is based on the AlphaFold predicted structures, crystal structures, cryo-EM, SAXS, connectivity and stereochemistry. The modeling was performed using [IMP](https://integrativemodeling.org) (*Integrative Modeling Platform*).

## Directory structure
1. [inputs](inputs/) : Contains the subdirectories for the input data used for the modeling.
2. [scripts](scripts/) : Contains all the scripts used for modeling, analysis and validation of the models.
3. [results](results/) : Contains the dominant cluster center model and the localization probability densities of the top cluster of the subcomplexes.

## Protocol
### Sampling
To perform sampling:
```
./run_modeling.sh <output_dir> <number_of_runs> <number_of_steps>
```

where, \
`output_dir` Path to the output directory \
`number_of_runs` Total number of sampling runs to perform\
`number_of_steps` Total number of 

### Analysis

Analysis was performed with the scripts in [IMP_Toolbox](https://github.com/isblab/IMP_Toolbox/tree/main/IMP_Toolbox/analysis)

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
* `aligned_*.pdb` : Atomic structures fitted to the cluster center model.

## Information

**Author(s):** Nisha Nandhini Shankar, Omkar Golatkar, Shruthi Viswanath, U. Venkatasubramanian1, Ragothaman M. Yennamalli\
**Date:**  21st August 2026\
**Publications:** Integrative structural modeling reveals the global architecture and conformational heterogeneity of _Acetivibrio thermocellus_ cellulosome
