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
from schrodinger.job.jobcontrol import launch_job
from schrodinger.application.jaguar.output import JaguarOutput
from schrodinger.structutils import analyze


def chunk_list(input_list, chunk_size):
    """Yield successive chunks of size `chunk_size` from `input_list`."""
    for i in range(0, len(input_list), chunk_size):
        yield input_list[i:i + chunk_size]



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Process MAE and input files")

    #parser.add_argument('-m', '--mae_path', type=str, required=True, help='Path to .mae file')
    parser.add_argument('-i', '--in_path', type=str, required=True, help='Path to .in file')
    parser.add_argument('-c', '--num_cpus', type=int, default=16, help='Number of CPU processors to use (default: 16) NEEDS TO BE DIVISIBLE BY NUM_PROCESS')
    parser.add_argument('-p', '--num_processes', type=int, default=4, help='Number of concurrent processes to spawn (default: 4)')
    args = parser.parse_args()

    in_path = os.path.realpath(args.in_path)
    mae_path = in_path.replace(".in",".mae")
    num_processes=num_processes = args.num_processes
    num_cpus = args.num_cpus

    mae_fname=os.path.basename(mae_path)
    molid_list = textscrape.get_molids(in_path)
    print(f"molIDs:{molid_list}")
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
            
            residues = list(mol.residue)
            res_batches = list(chunk_list(residues, args.num_processes))  # Break residues into batches of size n_cpus

            for batch_num, res_batch in enumerate(res_batches):
                job_inputs = []  # Collect input files for this batch

                for res in res_batch:
                    copy_st = mae_st.copy() # Make copy of original structure 
                    atom_indices = analyze.evaluate_asl(copy_st,f"{res.getAsl()} AND ( sidechain )") #evaluate ASL of residue sidechain
                    for i in atom_indices: # Set partial charge of all atoms in residue to 0
                        atom = copy_st.atom[i]
                        copy_st.atom[i].partial_charge = 0
                        copy_st.atom[i].solvation_charge = 0
                        atom.color = (0, 255, 0) # color atoms green for debug

                    res_num_str = str(res.resnum).zfill(5)
                    mae_copy_path = os.path.join(r_dir, f'{res_num_str}.mae')
                    in_copy_path = mae_copy_path.replace(".mae", ".in")

                    #write resnum.mae and resnum.in
                    with open(in_path, 'r') as f:
                        content = f.read().replace(os.path.basename(mae_path),os.path.basename(mae_copy_path))
                    with open(in_copy_path, 'w') as f:
                        f.write(content)
                    with structure.StructureWriter(mae_copy_path) as writer:
                        writer.append(copy_st)

                    job_inputs.append(os.path.basename(in_copy_path))  # Collect .in files for this batch
                
                os.chdir(r_dir)
                # Launch the batch as one job (can be modified to parallelize)
                run_command = ["jaguar", "run"] + job_inputs + ["-jobname", f"batch_{batch_num}", "-HOST", "localhost",'-PARALLEL',str(args.num_cpus)]
                print(run_command)
                job = launch_job(run_command)
                job.wait()
                os.chdir('..')
                for in_file in job_inputs:
                    try:
                        in_path = os.path.join(r_dir,in_file)
                        out_path = in_path.replace(".in", ".out")
                        result = output.QSiteOutput(out_path)
                        wavelength_nm = textscrape.extract_first_wavelength(out_path)
                        print(f"Total_E: {result.energy}")
                        print(f"Wavelength: {wavelength_nm}")
                        with open(os.path.join(p_dir,f"{res_num_str}.txt"),"w") as f:
                            f.write(f"{res_num_str}\t{result.energy}\t{wavelength_nm}")
                    except:
                        print(f'error reading {out_path}')

