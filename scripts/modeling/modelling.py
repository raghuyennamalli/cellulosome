"""
#############################################
##  cellulosome IMP modelling Script       ##
##                                         ##
#############################################
#
# Modeling script combining EM, distance restraint, proximity restraint
# 


"""
import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container
import RMF
import IMP.rmf

from IMP.pmi.restraints.basic import DistanceRestraint
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
import IMP.pmi.restraints.basic
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.macros

import IMP.pmi.topology
import ihm.cross_linkers
import IMP.display

import os
import sys
import time

from parameters import em_data_list_cellulosome
from parameters import distance_restraint_list_cellulosome

#Prints the selected random number
rng = IMP.random_number_generator()
print(rng)   # <class 'int'>
'''
#Random number can be fixed for reproducibility
seed_str = os.environ.get("SEED")
if seed_str is None:
    raise RuntimeError("SEED environment variable not set")
seed = int(seed_str)
IMP.random_number_generator.seed(seed)
print (f"random number set is {seed}")
'''
start_time = time.time()

#---------------------------
# 1. Define Input Data and Output Directories
#---------------------------
datadirectory = "../data/"
topology_file = datadirectory+"topology.txt"
output_directory = sys.argv[1]
output_index = sys.argv[2]

#--------------------------
# 2. Scoring Parameters
#--------------------------

#--------------
# ----- Sterochemistry and Physical Restraints
ev_weight = 1.0           # Weight of excluded volume restraint
connectivity_scale = 2.0     # weight of Connectivity restraint 2.0

#--------------------
# 3. Sampling Parameters
#--------------------
num_frames = int(sys.argv[3])  # Number of frames in MC run
num_mc_steps = 10         # Number of MC steps per frame
mc_temperature = 1.0      # Temperature for MC

# --- Simulated Annealing (sa)
#  - Alternates between two MC temperatures
sim_annealing = False    # If true, run simulated annealing
sa_min_temp_steps = 100 # Steps at min temp
sa_max_temp_steps = 20  # Steps at max temp
sa_temps = (1.0, 2.5)   # Sim annealing temperatures

# Replica Exchange (rex)
rex_temps = (1.0, 2.5)  # Temperature bounds for replica exchange

#------------------------------
# Here is where the real work begins...
#------------------------------

# Initialize model
m = IMP.Model()
outputobjects = [] # reporter objects (for stat files)

# Read in the topology file
# Specify the directory wheere the PDB files, fasta files and GMM files are
topology = IMP.pmi.topology.TopologyReader(topology_file,
                                  pdb_dir=datadirectory,
                                  fasta_dir=datadirectory,
                                  gmm_dir=datadirectory)

# Use the BuildSystem macro to build states from the topology file
bs = IMP.pmi.macros.BuildSystem(m)
 
# Each state can be specified by a topology file.
bs.add_state(topology)

# Build the system representation and degrees of freedom
root_hier, dof = bs.execute_macro(max_rb_trans=2.0, # 4.0 last worked 2.0
                                  max_rb_rot=0.3,   #0.3
                                  max_bead_trans=5.0, #5.0
                                  max_srb_trans=4.0,
                                  max_srb_rot=0.3)


root_hierarchy = root_hier

# Get copies, which are the immediate children of the root hierarchy
copy_particles = [child for child in root_hierarchy.get_children()]

# Optional: to print the entire hierarchy tree
def print_hierarchy(particle, level=0):
    print("  " * level + particle.get_name())
    for child in particle.get_children():
        print_hierarchy(child, level + 1)

#print_hierarchy(root_hierarchy)

def print_fragments(particle, level=0):
    name = particle.get_name()

    indent = "  " * level

    # check if particle is a fragment
    if IMP.atom.Fragment.get_is_setup(particle):
        frag = IMP.atom.Fragment(particle)
        residues = frag.get_residue_indexes()
        print(f"{indent}{name}  | residues: {list(residues)} | beads: 1")

    else:
        print(indent + name)

    for child in particle.get_children():
        print_fragments(child, level + 1)


#print_fragments(root_hierarchy)

#Selecting particles which are fitted to em map and temporarily disabling them to prevent their movement during shuffling
print("Selecting EM regions:")
em_particles = []

#cela.2, cela.3, cela.4: EM-fitted domains
for cp in [2,3,4]:
    #Enzymatic component
    sel = IMP.atom.Selection(root_hier, molecule=f"cela", 
                            copy_index=cp, residue_indexes=range(33,395))
    particles = sel.get_selected_particles()
    em_particles += particles
    print(f"  cela.{cp} [33-394]: {len(particles)} particles")
    
    #Dockerin  
    sel = IMP.atom.Selection(root_hier, molecule=f"cela", 
                            copy_index=cp, residue_indexes=range(413,477))
    particles = sel.get_selected_particles()
    em_particles += particles
    print(f"  cela.{cp} [413-476]: {len(particles)} particles")

# Cipa.0: EM-fitted domains 
cipa_ranges = [
    (560,704),   # coh3
    (724,867),   # coh4 
    (889,1032)   # coh5 
]
for start, end in cipa_ranges:
    sel = IMP.atom.Selection(root_hier, molecule="cipa", 
                            copy_index=0, residue_indexes=range(start, end))
    particles = sel.get_selected_particles()
    em_particles += particles 
    print(f"  cipa.0 [{start}-{end-1}]: {len(particles)} particles")

print(f"\nTotal EM particles: {len(em_particles)}")


# Connectivity keeps things connected along the backbone (ignores if inside same rigid body)

mols = IMP.pmi.tools.get_molecules(root_hier)
for mol in mols:
    molname=mol.get_name()
    IMP.pmi.tools.display_bonds(mol)
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol,scale=connectivity_scale, label=molname)
    cr.add_to_model()
    outputobjects.append(cr)

temp_hier = IMP.pmi.tools.input_adaptor(
    em_particles,
    pmi_resolution="all",
    flatten=True,
)

fixed_em_rigid_body, fixed_em_flexible_beads =  IMP.pmi.tools.get_rbs_and_beads(temp_hier)


# then restraints and continue sampling
# Create shuffle_frames directory inside the output directory
full_output_dir = f"{output_directory}{output_index}"
shuffle_dir = os.path.join(full_output_dir, "shuffle_frames")
os.makedirs(shuffle_dir, exist_ok=True)

# Construct the RMF filename
file_name_0 = f"initial_{output_directory}{output_index}.rmf3"
file_path_0 = os.path.join(shuffle_dir, file_name_0)

# Create and write the RMF file
fname_0 = RMF.create_rmf_file(file_path_0)   # create RMF file
IMP.rmf.add_hierarchy(fname_0, root_hier)  # add system hierarchy
IMP.rmf.save_frame(fname_0)                # write one frame
fname_0.close()                            # close file

#mild shuffling of em fitted particles
IMP.pmi.tools.shuffle_configuration(temp_hier,
                                    max_translation=10, 
                                    verbose=False,
                                    cutoff=4.0, 
                                    niterations=500)

# Construct the RMF filename
file_name_1 = f"mild_shuffle_{output_directory}{output_index}.rmf3"
file_path_1 = os.path.join(shuffle_dir, file_name_1)

# Create and write the RMF file
fname_1 = RMF.create_rmf_file(file_path_1)   # create RMF file
IMP.rmf.add_hierarchy(fname_1, root_hier)  # add system hierarchy
IMP.rmf.save_frame(fname_1)                # write one frame
fname_1.close()                            # close file

#shuffling of non-em particles
IMP.pmi.tools.shuffle_configuration(root_hier,
                                    excluded_rigid_bodies=fixed_em_rigid_body,
                                    max_translation=40, #40 
                                    verbose=False,
                                    cutoff=4.0, 
                                    niterations=500)
dof.optimize_flexible_beads(1000)

# Construct the RMF filename
file_name_2 = f"final_shuffle_{output_directory}{output_index}.rmf3"
file_path_2 = os.path.join(shuffle_dir, file_name_2)

# Create and write the RMF file
fname_2 = RMF.create_rmf_file(file_path_2)   # create RMF file
IMP.rmf.add_hierarchy(fname_2, root_hier)  # add system hierarchy
IMP.rmf.save_frame(fname_2)                # write one frame
fname_2.close()                            # close file

#-----------------------------------
# Define Scoring Function Components
#-----------------------------------

# Here we are defining a number of restraints on our system.
#  For all of them we call add_to_model() so they are incorporated into scoring
#  We also add them to the outputobjects list, so they are reported in stat files

# Excluded Volume Restraint
# To speed up this expensive restraint, we evaluate it at resolution 10
ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                         included_objects=root_hier,
                                         resolution=10)
ev.set_weight(ev_weight)
ev.add_to_model()
outputobjects.append(ev)

print ( "Excluded volume restraint appplied")
#  Electron Microscopy Restraint
#  The GaussianEMRestraint uses a density overlap function to compare model to data
#   First the EM map is approximated with a Gaussian Mixture Model (done separately)
#   Second, the components of the model are represented with Gaussians (forming the model GMM)
#   Other options: scale_to_target_mass ensures the total mass of model and map are identical
#                  slope: nudge model closer to map when far away
#                  weight: experimental, needed becaues the EM restraint is quasi-Bayesian
# Assume you have a dict mapping EM map filenames to the list of component names assigned:
print ("Starting EM particles selection loop...")

for em_params in em_data_list_cellulosome:

    density_particles = None

    for mol_name, res_range in em_params["molecules"].items():
        copy_idx = em_params.get("copy_idxs", {}).get(mol_name, 0)

        print(
            f"Mapping to component: {mol_name}, " 
            f"residue range: {res_range}, copy: {copy_idx}"
        )

        selection = IMP.atom.Selection(
            root_hier,
            molecule=mol_name,
            residue_indexes=list(res_range),
            representation_type=IMP.atom.DENSITIES,
            copy_index=copy_idx
        )

        if density_particles is None:
            density_particles = selection
        else:
            density_particles |= selection

    #MUST be inside this loop
    em_restraint = IMP.pmi.restraints.em.GaussianEMRestraint(
        densities=density_particles.get_selected_particles(),
        target_fn=em_params["gmm_file"],
        slope=em_params["slope"],
        scale_target_to_mass=True,
        weight=em_params["weight"]
    )

    em_restraint.set_label(em_params["label"])
    em_restraint.add_to_model()
    outputobjects.append(em_restraint)

print ("EM restraint applied")

#Distance restraint for cohesin 6, 7 and 8

print("Setting up distance restraints between cohesins 6,7 and 7,8..")

for distance_restraint_params in distance_restraint_list_cellulosome:
    distancemin = distance_restraint_params["target_distance"] - distance_restraint_params["threshold"]
    distancemax = distance_restraint_params["target_distance"] + distance_restraint_params["threshold"]

    tuple_selection1 = (distance_restraint_params["residue1_start"],
                        distance_restraint_params["residue1_end"],
                        distance_restraint_params["prot1"],
                        distance_restraint_params["copy1"])

    tuple_selection2 = (distance_restraint_params["residue2_start"],
                        distance_restraint_params["residue2_end"],
                        distance_restraint_params["prot2"],
                        distance_restraint_params["copy2"])
    print (f"Selected tuples: {tuple_selection1}, {tuple_selection2}")
    distance_restraint = DistanceRestraint(root_hier=root_hier,
                                           tuple_selection1=tuple_selection1,
                                           tuple_selection2=tuple_selection2,
                                           distancemin=distancemin,
                                           distancemax=distancemax,
                                           resolution=1.0,
                                           kappa=distance_restraint_params["kappa"],
                                           label=distance_restraint_params["label"],
                                           weight=distance_restraint_params["weight"])
    print(f"Creating DistanceRestraint between {distance_restraint_params['prot1']} residues {distance_restraint_params['residue1_start']}–{distance_restraint_params['residue1_end']} "
          f"and {distance_restraint_params['prot2']} residues {distance_restraint_params['residue2_start']}–{distance_restraint_params['residue2_end']}")
    print(f"Distance bounds: min={distancemin} Å, max={distancemax} Å; kappa={distance_restraint_params['kappa']}; label={distance_restraint_params['label']}")

    distance_restraint.add_to_model()
    outputobjects.append(distance_restraint)
    
print ("Distance restraint applied!")

#--------------------------
# Monte-Carlo Sampling
#--------------------------

if '--test' in sys.argv: num_frames=2000

# This object defines all components to be sampled as well as the sampling protocol

mc1=IMP.pmi.macros.ReplicaExchange(m,
                                   root_hier=root_hier,
                                   monte_carlo_sample_objects=dof.get_movers(),
                                   output_objects=outputobjects,
                                   monte_carlo_temperature=1.0,
#                                   simulated_annealing=False,
#                                   simulated_annealing_minimum_temperature=min(sa_temps),
#                                   simulated_annealing_maximum_temperature=max(sa_temps),
#                                   simulated_annealing_minimum_temperature_nframes=sa_min_temp_steps,
#                                   simulated_annealing_maximum_temperature_nframes=sa_max_temp_steps,
                                   replica_exchange_minimum_temperature=min(rex_temps),
                                   replica_exchange_maximum_temperature=max(rex_temps),
                                   number_of_best_scoring_models=0,
                                   monte_carlo_steps=num_mc_steps,
                                   number_of_frames=num_frames,
                                   global_output_directory=output_directory+str(output_index))

# Start Sampling
mc1.execute_macro()

print("started at: ", time.ctime(start_time))
print("ended at: ", time.ctime())
print('TOTAL TIME: ', round(time.time() - start_time, 2))
