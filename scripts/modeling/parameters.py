em_data_list_cellulosome = [
    {
        "name": "chimera_fitted_cohesin3_dockerin.pdb",
        "gmm_file": "../data/gmm/cohesin3_dockerin_15.txt",
        "slope": 0.000001, #0.000001
        "weight": 40.0,
        "label": "Cohesin_3",
        "molecules": {"cipa": range(560,704), 
                      "cela": range(413,478)},
        "copy_idxs": {"cipa":0, 
                      "cela": 2}
    },

    {
        "name": "chimera_fitted_cohesin4_dockerin.pdb",
        "gmm_file": "../data/gmm/cohesin4_dockerin_25.txt",
        "slope": 0.000001,
        "weight": 40.0,
        "label": "Cohesin_4",
        "molecules": {"cipa": range(719,867), 
                      "cela": range(413,478)},
        "copy_idxs": {"cipa": 0, 
                      "cela": 3}
    },

    {
        "name": "chimera_fitted_cohesin5_dockerin.pdb",
        "gmm_file": "../data/gmm/cohesin5_dockerin_55.txt",
        "slope": 0.000001,
        "weight": 40.0,
        "label": "Cohesin_5",
        "molecules": {"cipa": range(884,1032), 
                      "cela": range(413,478)},
        "copy_idxs": {"cipa":0, 
                      "cela": 4}
    },

    {
        "name": "chimera_fitted_1KWF_3",
        "gmm_file": "../data/gmm/enzyme3_5.txt",
        "slope": 0.000001,
        "weight": 40.0,
        "label": "1KWF_1",
        "molecules": {"cela": range(33,396)},
        "copy_idxs": {"cela": 2}
    },

    {
        "name": "chimera_fitted_1KWF_4",
        "gmm_file": "../data/gmm/enzyme4_5.txt",
        "slope": 0.000001,
        "weight": 40.0,
        "label": "1KWF_2",
        "molecules": {"cela": range(33,396)},
        "copy_idxs": {"cela": 3}
    },

    {
        "name": "chimera_fitted_1KWF_5",
        "gmm_file": "../data/gmm/enzyme5_5.txt",
        "slope": 0.000001,
        "weight": 40.0,
        "label": "1KWF_3",
        "molecules": {"cela": range(33,396)},
        "copy_idxs": {"cela": 4}
    }
]


distance_restraint_list_cellulosome = [
   {
    "prot1": "cipa",
    "copy1": 0,
    "residue1_start": 967,
    "residue1_end": 967,
    "prot2": "cipa",
    "copy2": 0,
    "residue2_start": 1132,
    "residue2_end": 1132,
    "target_distance": 58.743,  #distance found from chimerax
    "threshold": 5.0,           #upper and lower bound to be added to the target distance
    "kappa": 0.5, #1.0
    "weight": 1.0, #1.0
    "label": "coh5_0_967_coh6_0_1132_COM_Distance"
},
   {
    "prot1": "cipa",
    "copy1": 0,
    "residue1_start": 1132,
    "residue1_end": 1132,
    "prot2": "cipa",
    "copy2": 0,
    "residue2_start": 1297,
    "residue2_end": 1297,
    "target_distance": 58.743,
    "threshold": 5.0,           #upper and lower bound to be added to the target distance
    "kappa": 0.5, #1.0
    "weight": 1.0, #1.0
    "label": "coh6_0_1132_coh7_0_1297_COM_Distance"
},
    {
    "prot1": "cipa",
    "copy1": 0,
    "residue1_start": 1297,
    "residue1_end": 1297,
    "prot2": "cipa",
    "copy2": 0,
    "residue2_start": 1462,
    "residue2_end": 1462,
    "target_distance": 58.743,
    "threshold": 5.0,
    "kappa": 0.5,
    "weight": 1.0,
    "label": "coh7_0_1297_coh8_0_1462_COM_Distance"
},

   {
    "prot1": "cipa",
    "copy1": 0,
    "residue1_start": 1462,
    "residue1_end": 1462,
    "prot2": "cipa",
    "copy2": 0,
    "residue2_start": 1626,
    "residue2_end": 1626,
    "target_distance": 58.743,
    "threshold": 5.0,
    "kappa": 0.5,
    "weight": 1.0,
    "label": "coh8_0_1462_coh9_0_1626_COM_Distance"
},

    {
    "prot1": "cipa",
    "copy1": 0,
    "residue1_start": 29,
    "residue1_end": 29,
    "prot2": "cipa",
    "copy2": 0,
    "residue2_start": 1853,
    "residue2_end": 1853,
    "target_distance": 575,    # D_max from SAXS
    "threshold": 20.0,           # Threshold ±20 Å
    "kappa": 0.5,
    "weight": 1,
    "label": "cellulosome_end_to_end_Distance"
}
]
