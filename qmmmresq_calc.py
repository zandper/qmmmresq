#!/usr/bin/env python3
from schrodinger import structure
from schrodinger.application.qsite import output
#from schrodinger.test import mmshare_data_fil

import os
import re
import shutil
import subprocess
import multiprocessing
import argparse
from utils import textscrape

def calc_single_point_residue(res, mol, mae_st, r_dir, mae_fname, in_path, p_dir):
    copy_st = mae_st.copy() # Make copy of original structure 
    atom_list = res._atoms # Get index of all atoms in residue
    for i in atom_list: # Set partial charge of all atoms in residue to 0 
        copy_st.atom[i].partial_charge = 0
        copy_st.atom[i].solvation_charge = 0 #ASK ABOUT THE SOLVATION CHARGE

    res_num_str = str(res.resnum).zfill(5)
    print(res_num_str)
    res_dir=os.path.join(r_dir,res_num_str)
    os.makedirs(res_dir,exist_ok=True)
    mae_copy_path = os.path.join(res_dir,mae_fname)
    in_copy_path = mae_copy_path.replace(".mae",".in")
    shutil.copyfile(in_path, in_copy_path)
    with structure.StructureWriter(mae_copy_path) as writer: #write modified .mae file
        writer.append(copy_st)
    cmd_dir = os.path.realpath(res_dir)
    p = subprocess.Popen(['qsite', '-WAIT', os.path.basename(in_copy_path)], cwd=cmd_dir) #run qsite calculation
    p.wait()
    out_path = mae_copy_path.replace(".mae",".out")
    result = output.QSiteOutput(out_path)
    wavelength_nm = textscrape.extract_first_wavelength(out_path)
    print(f"Total_E: {result.energy}")
    print(f"Wavelength: {wavelength_nm}")


    with open(os.path.join(p_dir,f"{res_num_str}.txt"),"w") as f:
        f.write(f"{res_num_str}\t{result.energy}\t{wavelength_nm}")
    #shutil.rmtree(res_dir) 


def process_all_residues(mol, mae_st, r_dir, mae_fname, in_path, p_dir):
    # Use multiprocessing Pool to parallelize the work
    with multiprocessing.Pool(processes=num_processes) as pool:
        # Pass all arguments required for each residue to the pool
        args = [(res, mol, mae_st, r_dir, mae_fname, in_path, p_dir) for res in mol.residue]
        pool.starmap(calc_single_point_residue, args)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Process MAE and input files")

    #parser.add_argument('-m', '--mae_path', type=str, required=True, help='Path to .mae file')
    parser.add_argument('-i', '--in_path', type=str, required=True, help='Path to .in file')
    parser.add_argument('-p', '--num_processes', type=int, default=12, help='Number of CPU processes to use (default: 12)')
    args = parser.parse_args()

    in_path = os.path.realpath(args.in_path)
    mae_path = in_path.replace(".in",".mae")
    num_processes=12

    mae_fname=os.path.basename(mae_path)
    
    molid_list = textscrape.get_molids(in_path)

    print(molid_list)

    mae_st = structure.StructureReader.read(mae_path)

    p_dir="parameters"
    os.makedirs(p_dir, exist_ok=True)
    r_dir="residues"
    os.makedirs(r_dir, exist_ok=True)


    #iterate through all molecules in structure
    for mol in mae_st.molecule:
        #Skip molecules which are in qmm region
        if mol.number in molid_list:
            print(f"molecule #{mol.number} in qmm region, skipping molecule")
        else:
            #Asynchronously calculate for all residues in a molecule
            process_all_residues(mol, mae_st, r_dir, mae_fname, in_path, p_dir)

        """
        for res in mol.residue:
            calc_single_point_residue(res, mol, mae_st, r_dir, mae_fname, in_path, p_dir)
            exit()
        """
            #DO THIS PART ASYNCHRONOUSLY

        #t = structure._AtomCollection(mae_st, atom_list)
        """
        print(t[0])
        exit()
        st_temp = mae_st.copy()
        st_temp.property = st_temp.property
        print(st_temp.property)
        atom_list= res._atoms
        print(mae_st.atom[atom_list])
        """
        #nres = build.neutralize_structure(res)
        #print(dir(atoms))
        #print(res.getCode())
        #print(res.resnum) #resnum is the mol
        #atom_list= res._atoms
        #set atom charges to 0
        #print(atom_list)
        #mae_st()
        #set charge of residue equal to 0
        
        #save .mae file
        #save .in file (needs to have matching name)
        #run calc

        #print(dir(res))

        #temp_mae_st = mae_st.neutralize residue

            

#USE THIS FOR MEASURING PROXIMITY TO CHROMOPHORE
#schrodinger.structutils.measure.get_shortest_distance(st, atoms=None, st2=None, cutoff=inf)


"""
residues = residues not in qm region x
for i in residues: x
    set residue to neutral charge x
    calc single point energy
    write energy
    chromophore distance
"""

"""
plot delta energy vs residue (sort by chromophore distance)
"""