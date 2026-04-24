import re
import os
import sys
import time
import IMP
import IMP.rmf
import IMP.atom
import RMF
import multiprocessing as mp
import numpy as np
import tqdm
import pickle
from collections import defaultdict

def parse_rmf_into_radii_and_coords(location, burn_in):
    # Extracts coordinates and radii from an RMF file
    m = IMP.Model()
    rmf_file = RMF.open_rmf_file_read_only(location)
    n = rmf_file.get_number_of_frames()
    print(f'Number of Models in the file {location} is {n} with the initial {burn_in} removed')
    h = IMP.rmf.create_hierarchies(rmf_file, m)[0]
    model_wise_coords = []
    molecule_wise_coords = []
    # to be converted later to a dictionary with molecule name -> models x particle coordinates
    particle_radii = []  # all particle radii
    molecule_wise_radii = defaultdict(list)  # molecule name -> radius of all the particles in order
    molecule_wise_names = defaultdict(list)  # molecule name -> names of all the particles in order

    for i in tqdm.trange(burn_in, n):
        IMP.rmf.load_frame(rmf_file, i)
        m.update()
        all_particles = [x.get_particle() for x in IMP.atom.get_leaves(h)]
        all_coords = []
        for p in all_particles:
            xyzr = IMP.core.XYZR(p)
            all_coords.append(np.array(xyzr.get_coordinates()))
            if i == burn_in:
                particle_radii.append(xyzr.get_radius())
        molecules = h.get_children()[0].get_children()
        molecule_names = [x.get_name() for x in molecules]
        molecule_coords = defaultdict(list)
        
        for name, mol in zip(molecule_names, molecules):
            sub_array = []
            radii = []
            for x in IMP.atom.get_leaves(mol):
                xyzr = IMP.core.XYZR(x.get_particle())
                sub_array.append(np.array(xyzr.get_coordinates()))
                radii.append(xyzr.get_radius())
            
            molecule_coords[name] += sub_array
            if i == burn_in:
                molecule_wise_radii[name] += radii[:]
                molecule_wise_names[name] += [x.get_name() for x in IMP.atom.get_leaves(mol)]
        
        model_wise_coords.append(all_coords)
        molecule_wise_coords.append(molecule_coords)
    model_wise_coords = np.array(model_wise_coords)
    molecule_wise_coords_dict = dict()
    for i in molecule_wise_coords[0]:
        molecule_wise_coords_dict[i] = np.array([j[i] for j in molecule_wise_coords])
    print(f'Finished loading data from {location}')
    return model_wise_coords, molecule_wise_coords_dict, molecule_wise_names, molecule_wise_radii, particle_radii

def sort_only_columns(matrix):
    # to keep the columns intact and sort them
    return np.array(sorted(matrix.T, key=lambda x: x.tolist())).T
    
def foo_rmf_parser(r, ind, burn_in, q_save):
    # Worker function for parallel RMF parsing
    m1, m2, names, p1, p2 = parse_rmf_into_radii_and_coords(r, burn_in)
    q_save.put((m2, p1, ind))
    return names

def foo_saver(total, q, save_path):
    m_aggregate = None
    p_aggregate = None
    for i in range(total):
        m, p, name = q.get()
        if p_aggregate is None:
            p_aggregate = p
        else:
            for k in p:
                n = len(p[k])
                assert all([p[k][j] == p_aggregate[k][j] for j in range(n)]), f'Radii are inconsistent: {i}'
        if m_aggregate is None:
            m_aggregate = m
        else:
            for k in m_aggregate:
                m_aggregate[k] = np.vstack([m_aggregate[k], m[k]])
    with open(f'{save_path}/saved_data', 'wb') as f:
        pickle.dump([m_aggregate, p_aggregate], f)
    print('Data saved')
    
#checked-----------------------------
def process_names(names):
    # Maps bead/fragment names to residue indices
    indices = defaultdict(list)
    residues = defaultdict(list)
    bead_sort = defaultdict(list)
    reg = '[^0-9]*([0-9]+)[^0-9]+([0-9]+)[^0-9]*'
    f = open('running_temp_name_fragment_map.txt', 'w')
    for k in names:
        already_done = set()
        f.write(f'{k}\n')
        for i, name in enumerate(names[k]):
            if name in already_done: 
                continue
            already_done.add(name)
            if 'bead' in name:
                #modified
                match = re.search(reg, name)
                assert match, f'Regex pattern failed to match {name}'
                n1 = int(match.group(1))
                n2 = int(match.group(2))
                residues[k] += list(range(n1, n2 + 1))
                indices[k] += [i for _ in range(n1, n2 + 1)]
                f.write(f'{n1} -> {n2}\n')
                bead_sort[k].append(n1)
            elif 'Fragment' in name:
                match = re.search(reg, name)
                assert match, f'Regex pattern failed to match {name}'
                n1 = int(match.group(1))
                n2 = int(match.group(2))
                residues[k] += list(range(n1, n2))
                indices[k] += [i for _ in range(n1, n2)]
                f.write(f'{n1} -> {n2 - 1}\n')
                bead_sort[k].append(n1)
            #modified
            elif name.strip().isdigit():
                # single residue bead like '33'
                res_num = int(name.strip())
                residues[k] += [res_num]
                indices[k] += [i]
                f.write(f'{res_num}\n')
                bead_sort[k].append(res_num)
            else:
                assert False, f'Neither bead nor fragment: {name} in {k}'
    f.close()
    return indices, residues, bead_sort

def get_all_data_xyzr(list_of_rmfs, npool=None, save_path='.', burn_in=0):
    if npool is None:
        npool = min(len(list_of_rmfs), os.cpu_count())
    with mp.Pool(npool) as p:
        print('Setting up.')
        manager = mp.Manager()
        q3 = manager.Queue()
        cm_save = f'{save_path}/contact_maps'
        data_save = f'{save_path}/extracted_xyzr'
        if not os.path.isdir(data_save):
            os.mkdir(data_save)
        if not os.path.isdir(cm_save):
            os.mkdir(cm_save)
            args = [(list_of_rmfs[r], r, burn_in, q3) for r in range(len(list_of_rmfs))]
            p3 = mp.Process(target=foo_saver, args=(len(list_of_rmfs), q3, data_save), daemon=True)
            p3.start()
            print('Saving molecule-wise coordinates and radii.')
            temp = p.starmap(foo_rmf_parser, args)
            assert all(temp), 'Some of the foo_workers did not end properly while calculating parsing rmfs'
            names = temp[0]
            print('Finished!')
            indices, residues, bead_sort = process_names(names)
            del temp
            print('Waiting for processes to die')
            wait_time_start = time.time()
            while (time.time() - wait_time_start) < 60:
                p3.join(5)
                if not p3.is_alive():
                    print('Processes have ended properly.')
                    break
            if p3.is_alive():
                print('Processes have not ended. Exiting the function.')
                print(f'\tProcess p3: {p3.is_alive()}')
            with open(f'{cm_save}/names_indices', 'wb') as f:
                pickle.dump([indices, residues, bead_sort, names], f)
        else:
            with open(f'{cm_save}/names_indices', 'rb') as f:
                indices, residues, bead_sort, names = pickle.load(f)
            
    print('XYZR extraction and metadata generation complete.')

if __name__ == '__main__':
    # usage: python script.py <input_dir> <unused> <save_path> <threads> <burn_in> <cluster_id>
    input_dir = sys.argv[1]
    save_path = sys.argv[3]
    n_pool = int(sys.argv[4])
    burn_in = int(sys.argv[5])
    cluster_number = int(sys.argv[6])
    
    location = [x for x in os.listdir(sys.argv[1]) if re.search(f'[AB]_gsm_clust{cluster_number}.rmf3', x)]
    rmfs = [f'{sys.argv[1]}/{x}' for x in location]
    get_all_data_xyzr(rmfs, n_pool, save_path, burn_in)
