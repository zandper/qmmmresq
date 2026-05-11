#!/usr/bin/env python3

#NEED TO ADD USE_MAE_CHARGES TO .IN FILE
from schrodinger import structure
from schrodinger.structutils import analyze, measure
from schrodinger.application.qsite import output
from schrodinger.application.qsite.input import QSiteInput
#from schrodinger.test import mmshare_data_fil
import os
import re
import shutil
import subprocess
import multiprocessing
import argparse
import glob
from utils import textscrape


def get_nearby_mol_res(mae_st, qm, cutoff_angstrom):
    """
    Returns a sorted list of (molnum, resnum) for residues
    with at least one atom within cutoff_angstrom of ligand ASL.
    Excludes atoms in the ligand itself.
    """
    # Ligand atoms
    qm_atoms = set(analyze.evaluate_asl(mae_st, qm))
    if not qm_atoms:
        raise ValueError(f"No atoms selected by ASL: '{qm}'")

    # All atoms within cutoff of ligand
    nearby_atoms = set(measure.get_atoms_close_to_subset(mae_st, qm_atoms, cutoff_angstrom))
    nearby_atoms -= qm_atoms     # Exclude QM region
    # Map atoms to residues
    residues = set()
    for idx in nearby_atoms:
        atom = mae_st.atom[idx]
        residues.add((atom.molecule_number, atom.resnum))

    return sorted(residues)


def calc_single_point_residue(molnum, resnum, mae_path, r_dir, p_dir, n_cpu):

    # ---- Build file identidfiers----
    mae_fname = os.path.basename(mae_path)
    frame_num = mae_fname.split("_")[0]
    res_num_str = f"{str(frame_num).zfill(6)}_{str(molnum).zfill(6)}_{str(resnum).zfill(6)}"
    res_dir = os.path.join(r_dir, res_num_str)
    in_path = mae_path.replace("_geopt.01.mae","_spe.in")
    mae_copy_path = os.path.join(res_dir, mae_fname)
    out_path = os.path.join(res_dir,os.path.basename(in_path).replace("01.in", "02.out"))

    # ---- Skip if completed ----
    if os.path.exists(out_path):
        try:
            with open(out_path) as f:
                text = f.read()
            if "completed on" in text.lower():
                print(f"[SKIP] {res_num_str} already completed.")
                return
        except Exception:
            print(f"[WARNING] Could not validate existing .out, rerunning.")

    print(f"calculating {res_num_str}")
    os.makedirs(res_dir, exist_ok=True)

    # ---- Load structure fresh ----
    mae_st = structure.StructureReader.read(mae_path)

    # Select target residue
    target_mol = next(m for m in mae_st.molecule if m.number == molnum)
    target_res = next(r for r in target_mol.residue if r.resnum == resnum)

    sidechain_ASL = f"{target_res.getAsl()} AND NOT ( backbone )"
    atom_indices = analyze.evaluate_asl(mae_st, sidechain_ASL)

    for i in atom_indices:
        atom = mae_st.atom[i]
        atom.partial_charge = 0
        atom.solvation_charge = 0
        atom.color = (0, 255, 0)

    # ---- Prepare input files ----
    in_copy_path = mae_copy_path.replace(".mae", ".in")
    shutil.copyfile(in_path, in_copy_path)
    in_copy_path = os.path.realpath(in_copy_path)
    textscrape.add_mae_charges_yes(in_copy_path) # add use mae charges argument to infile

    with structure.StructureWriter(mae_copy_path) as writer:
        writer.append(mae_st)

    # ---- Run QSite ----
    cmd = ['qsite', '-WAIT', '-HOST', 'localhost','-PARALLEL', str(int(n_cpu)),os.path.basename(in_copy_path)]
    print(cmd)
    p = subprocess.Popen(cmd, cwd=res_dir)
    p.wait()

    # ---- Post-processing ----
    print(f"OUT_PATH:{out_path}")
    try:
        result = output.QSiteOutput(out_path)
        wavelength_nm = textscrape.extract_first_wavelength(out_path)

        print(f"Total_E: {result.energy}, Wavelength: {wavelength_nm}")

        summary_path = os.path.join(p_dir, f"{res_num_str}.txt")

        with open(summary_path, "w") as f:
            f.write(f"{molnum}\t{resnum}\t{result.energy}\t{wavelength_nm}\n")

    except Exception as e:
        print(f"Error reading {out_path}: {e}")


def process_all_residues(mae_path, r_dir, p_dir, num_processes, n_cpu, molnum_resnum_list):
    
    args = [(molnum, resnum, mae_path, r_dir, p_dir, n_cpu) for molnum, resnum in molnum_resnum_list]
    
    print(f"Selected residues: {len(args)}")
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.starmap(calc_single_point_residue, args)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Process MAE and input files")
    parser.add_argument('-m', '--mae_folder', type=str, required=True, help='Path to folder containing .mae files')
    parser.add_argument('-c', '--num_cpus', type=int, default=16, help='Number of CPU processors to use per process')
    parser.add_argument('-p', '--num_processes', type=int, default=4, help='Number of concurrent processes to spawn (default: 4)')
    parser.add_argument('-qa','--qasl',type=str,default=None,help='Manually define residues in QM')
    parser.add_argument('-d','--distance',type=float, default=None,help='Limit range from QM region')

    args = parser.parse_args()
    num_processes = args.num_processes
    n_cpu = args.num_cpus
    
    mae_paths = [f for f in glob.glob(os.path.join(args.mae_folder, "*_geopt.01.mae"))]
    first_st = structure.StructureReader.read(mae_paths[0])
    molnum_resnum_list = get_nearby_mol_res(first_st, args.qasl, args.distance)
    print(f"evaluating {len(molnum_resnum_list)} residues")

    p_dir="parameters"
    os.makedirs(p_dir, exist_ok=True)
    r_dir="residues"
    os.makedirs(r_dir, exist_ok=True)
    counter = 0
    for mae_path in mae_paths:
            process_all_residues(mae_path, r_dir, p_dir, num_processes, n_cpu,molnum_resnum_list)


    #iterate through all molecules in structure

#        #Skip molecules which are in qmm region
#        if mol.number in molid_list:
#            print(f"molecule #{mol.number} in qmm region, skipping molecule")
#        else:
#            exit()
            # Asynchronous calculation
            #process_all_residues(mol, mae_st, r_dir, mae_fname, in_path, p_dir, num_processes, n_cpu)

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