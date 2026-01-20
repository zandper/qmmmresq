#!/usr/bin/env python3

#NEED TO ADD USE_MAE_CHARGES TO .IN FILE
from schrodinger import structure
from schrodinger.structutils import analyze
from schrodinger.application.qsite import output
#from schrodinger.test import mmshare_data_fil

import os
import re
import shutil
import subprocess
import multiprocessing
import argparse
from utils import textscrape

def calc_single_point_residue(res, mol, mae_st, r_dir, mae_fname, in_path, p_dir, n_cpu):
    copy_st = mae_st.copy() # Make copy of original structure 
    sidechain_ASL = f"{res.getAsl()} AND ( sidechain )" #Generate ASL of residue sidechain
    atom_indices = analyze.evaluate_asl(copy_st,sidechain_ASL) #evaluate ASL
    for i in atom_indices: # Set partial charge of all atoms in residue to 0
        atom= copy_st.atom[i]
        #print(f"atom: {i}")
        #print(f'type: {atom.atom_name}')
        copy_st.atom[i].partial_charge = 0
        copy_st.atom[i].solvation_charge = 0
        atom.color = (0, 255, 0) # color atoms green for debug

    res_num_str = str(res.resnum).zfill(5)
    print(f"calculating {res_num_str}")
    res_dir=os.path.join(r_dir,res_num_str)
    os.makedirs(res_dir,exist_ok=True)
    mae_copy_path = os.path.join(res_dir,mae_fname)
    in_copy_path = mae_copy_path.replace(".mae",".in")
    shutil.copyfile(in_path, in_copy_path)
    with structure.StructureWriter(mae_copy_path) as writer: #write modified .mae file
        writer.append(copy_st)
    cmd_dir = os.path.realpath(res_dir)
    cmd =['qsite', '-WAIT','-HOST','localhost','-PARALLEL', str(int(n_cpu)), os.path.basename(in_copy_path)]
    #print(cmd)
    p = subprocess.Popen(cmd, cwd=cmd_dir) #run qsite calculation
    p.wait()
    try:
        out_path = mae_copy_path.replace(".mae",".out")
        result = output.QSiteOutput(out_path)
        print(result)
        wavelength_nm = textscrape.extract_first_wavelength(out_path)
        print(f"Total_E: {result.energy}")
        print(f"Wavelength: {wavelength_nm}")

        with open(os.path.join(p_dir,f"{res_num_str}.txt"),"w") as f:
            f.write(f"{res_num_str}\t{result.energy}\t{wavelength_nm}")
            #shutil.rmtree(res_dir)
    except Exception as e:
        print(f"Error reading {out_path}: {e}")


def process_all_residues(mol, mae_st, r_dir, mae_fname, in_path, p_dir,num_processes,n_cpu):
    # Use multiprocessing Pool to parallelize the work
    with multiprocessing.Pool(processes=num_processes) as pool:
        # Pass all arguments required for each residue to the pool
        args = [(res, mol, mae_st, r_dir, mae_fname, in_path, p_dir,n_cpu) for res in mol.residue]
        pool.starmap(calc_single_point_residue, args)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Process MAE and input files")
    parser.add_argument('-i', '--in_path', type=str, required=True, help='Path to .in file')
    parser.add_argument('-c', '--num_cpus', type=int, default=16, help='Number of CPU processors to use (default: 16) NEEDS TO BE DIVISIBLE BY NUM_PROCESS')
    parser.add_argument('-p', '--num_processes', type=int, default=4, help='Number of concurrent processes to spawn (default: 4)')
    args = parser.parse_args()

    in_path = os.path.realpath(args.in_path)
    textscrape.add_mae_charges_yes(in_path)
    mae_path = in_path.replace(".in",".mae")
    mae_fname=os.path.basename(mae_path)
    molid_list = textscrape.get_molids(in_path)
    #print(molid_list)
    num_processes = args.num_processes
    n_cpu = args.num_cpus/num_processes
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

            # Asynchronous calculation
            process_all_residues(mol, mae_st, r_dir, mae_fname, in_path, p_dir, num_processes, n_cpu)

            #Synchronous calculation KEEP FOR DEBUG
            #for res in mol.residue:
            #    print(res.resnum)
            #    calc_single_point_residue(res, mol, mae_st, r_dir, mae_fname, in_path, p_dir,n_cpu)
            #    exit()
            #    pass





            

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