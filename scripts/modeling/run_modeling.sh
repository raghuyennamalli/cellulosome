#!/bin/bash

# Usage:
# ./run_modeling.sh out_dir n n_steps

py=/home/ragothaman/miniconda3/envs/IMP/bin/python

out_dir=$1
n=$2
n_steps=$3

# List of seeds for each replica
#seeds=(28 5 23 12 17 33 46 57 68 79 96 132 124 1255 151 515 124 112 563 234 236 982 232 1231 14124 1413241 2343534 5435453 64565634 32 342 325352)  

# Check that there are enough seeds
#if [ ${#seeds[@]} -lt $n ]; then
#    echo "Error: Not enough seeds for $n replicas. Add more seeds to the list."
#    exit 1
#fi

for x in $( seq 1 $n )
do
    # Get the seed for this replica from the list (0-indexed)
    #seed=${seeds[$((x-1))]}
    #echo "Running job $x with seed $seed..."
    echo "Running job $x..." 

    # Running without replica exchange
    #$py modelling.py $out_dir $x $n_steps > ./output_files/${out_dir}$x.out 2>&1

    # Running with replica exchange
    #SEED="$seed" mpirun -np 4 $py modelling_random_number_set.py $out_dir $x $n_steps > ./output_files/${out_dir}$x.out 2>&1
    mpirun -np 6 $py modelling.py $out_dir $x $n_steps > ./output_files/${out_dir}$x.out 2>&1

    #echo "Job $x with seed $seed completed."
    echo "Job $x completed."
done



