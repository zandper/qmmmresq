#!/usr/bin/env python3

#NEED TO ADD USE_MAE_CHARGES TO .IN FILE
from schrodinger import structure
from schrodinger.structutils import analyze, measure
from schrodinger.application.qsite import output
from schrodinger.application.qsite.input import QSiteInput
#from schrodinger.test import mmshare_data_fil
import os
import re
import glob
import shutil
import subprocess
import multiprocessing
import argparse
import utils.textscrape
import utils.qm_connectivity


#REMEMBER TO USE (chain.name A) FOR PROTEIN ASL
#REMEMBER TO USE (chain.name A) FOR PROTEIN ASL
#REMEMBER TO USE (chain.name A) FOR PROTEIN ASL
#REMEMBER TO USE (chain.name A) FOR PROTEIN ASL
#REMEMBER TO USE (chain.name A) FOR PROTEIN ASL
#REMEMBER TO USE (chain.name A) FOR PROTEIN ASL
#REMEMBER TO USE (chain.name A) FOR PROTEIN ASL
#REMEMBER TO USE (chain.name A) FOR PROTEIN ASL
#REMEMBER TO USE (chain.name A) FOR PROTEIN ASL

#alternatively accept molid if hcaps not present
def get_nearby_mol_res(mae_st, qm_asl, cutoff_angstrom=None, manual_asl=None, protein_only=False):
    """Get nearby molecule/residue information around QM region."""
    qm_atoms = set(analyze.evaluate_asl(mae_st, qm_asl))
    if not qm_atoms:
        raise ValueError(f"No atoms selected by ASL: '{qm_asl}'")
    
    if manual_asl:
        selected_atoms = set(analyze.evaluate_asl(mae_st, manual_asl))
        if not selected_atoms:
            raise ValueError(f"No atoms selected by manual ASL: '{manual_asl}'")
        nearby_atoms = selected_atoms
    elif cutoff_angstrom is not None:
        nearby_atoms = set(measure.get_atoms_close_to_subset(mae_st, qm_atoms, cutoff_angstrom))
        nearby_atoms -= qm_atoms
    else:
        raise ValueError("Either cutoff_angstrom or manual_asl must be provided")
    
    if protein_only:
        solvent_atoms = set(analyze.evaluate_asl(mae_st, "solvent"))
        membrane_atoms = set(analyze.evaluate_asl(mae_st, "membrane"))
        salt_atoms = set(analyze.evaluate_asl(mae_st, '(ions) AND ((chain.name " "))'))
        nearby_atoms -= solvent_atoms
        nearby_atoms -= membrane_atoms
        nearby_atoms -= salt_atoms

    return sorted({
        (mae_st.atom[i].molecule_number, mae_st.atom[i].resnum)
        for i in nearby_atoms
    })


def prepare_residue_files(molnum, resnum, mae_path, in_path, target_dir):
    """
    Prepare .mae and .in files for a single residue.
    Returns tuple (res_num_str, in_file_path) or None if failed.
    """
    res_num_str = f"{str(molnum).zfill(6)}_{str(resnum).zfill(6)}"
    mae_copy_path = os.path.join(target_dir, f"{res_num_str}.mae")
    in_copy_path = os.path.join(target_dir, f"{res_num_str}.in")
    
    try:
        # Load structure fresh
        mae_st = structure.StructureReader.read(mae_path)
        
        # Select target residue
        target_mol = next(m for m in mae_st.molecule if m.number == molnum)
        target_res = next(r for r in target_mol.residue if r.resnum == resnum)
        
        # Neutralize sidechain
        sidechain_ASL = f"{target_res.getAsl()} AND NOT ( backbone )"
        atom_indices = analyze.evaluate_asl(mae_st, sidechain_ASL)
        
        for i in atom_indices:
            atom = mae_st.atom[i]
            atom.partial_charge = 0 #Set sidechain charges to zero
            atom.solvation_charge = 0 #Set sidechain charges to zero
            atom.color = (0, 255, 0) #Color target sidechain for debugging
        
        # Prepare input files
        with open(in_path, 'r') as f:
            in_content = f.read()
        in_content = re.sub(r'MAEFILE:\s*\S+', f'MAEFILE: {res_num_str}.mae', in_content)
        
        with open(in_copy_path, 'w') as f:
            f.write(in_content)
        
        with structure.StructureWriter(mae_copy_path) as writer:
            writer.append(mae_st)
        
        return (res_num_str, in_copy_path)
        
    except Exception as e:
        print(f"Failed to prepare {res_num_str}: {e}")
        return None


def process_result(out_path, molnum, resnum, p_dir, native_lambda):
    """Process a single output file and write summary."""
    res_num_str = f"{str(molnum).zfill(6)}_{str(resnum).zfill(6)}"
    summary_path = os.path.join(p_dir, f"{res_num_str}.txt")
    
    if os.path.exists(summary_path):
        return True
    
    try:
        result = output.QSiteOutput(out_path)
        wavelength_nm = float(utils.textscrape.extract_first_wavelength(out_path))
        res_contrib = native_lambda - wavelength_nm
        print(f"Residue {res_num_str}: Total_E: {result.energy}, Contribution: {res_contrib}, Wavelength: {wavelength_nm}")
        
        with open(summary_path, "w") as f:
            f.write(f"{molnum}\t{resnum}\t{result.energy}\t{res_contrib}\t{wavelength_nm}\n")
        return True
    except Exception as e:
        print(f"Failed to process {res_num_str}: {e}")
        return False


def calc_single_point_residue(molnum, resnum, mae_path, r_dir, in_path, p_dir, n_cpu, native_lambda, max_retries=3):
    """Calculate single point residue contribution using QSite directly."""
    res_num_str = f"{str(molnum).zfill(6)}_{str(resnum).zfill(6)}"
    summary_path = os.path.join(p_dir, f"{res_num_str}.txt")
    
    if os.path.exists(summary_path):
        print(f"[SKIP] {res_num_str} already completed.")
        return
    
    for attempt in range(1, max_retries + 1):
        print(f"calculating {res_num_str} (attempt {attempt}/{max_retries})")

        
        try:
            # Prepare files
            prepared = prepare_residue_files(molnum, resnum, mae_path, in_path, r_dir)
            if not prepared:
                raise Exception("Failed to prepare files")
            
            res_num_str, in_copy_path = prepared
            out_path = in_copy_path.replace(".in", ".out")
            
            # Run QSite
            cmd = ['qsite', '-WAIT', '-HOST', 'localhost', '-PARALLEL', str(int(n_cpu)), os.path.basename(in_copy_path)]
            print(f"Running in {r_dir}: {' '.join(cmd)}")
            p = subprocess.Popen(cmd, cwd=r_dir)
            p.wait()
            
            if p.returncode != 0:
                raise Exception(f"QSite failed with return code {p.returncode}")
            
            # Process result
            if process_result(out_path, molnum, resnum, p_dir, native_lambda):
                #shutil.rmtree(res_dir)
                return
            
        except Exception as e:
            print(f"Attempt {attempt} failed: {e}")
            if attempt == max_retries:
                print(f"Giving up on {res_num_str}")
            else:
                print(f"removing{r_dir}")
                #shutil.rmtree(r_dir, ignore_errors=True)


def process_all_residues(mae_path, r_dir, in_path, p_dir, num_processes, n_cpu, native_lambda, molnum_resnum_list):
    """Process all residues using Python multiprocessing."""
    args = [(molnum, resnum, mae_path, r_dir, in_path, p_dir, n_cpu, native_lambda) 
            for molnum, resnum in molnum_resnum_list]
    
    print(f"Selected residues: {len(args)}")
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.starmap(calc_single_point_residue, args)


def calc_jaguar_parallel(mae_path, r_dir, in_path, p_dir, n_cpu, native_lambda, molnum_resnum_list):
    """Run calculations using Jaguar's parallel mode."""
    print("Preparing input files for all residues...")
    
    # Prepare all jobs
    in_files = []
    for molnum, resnum in molnum_resnum_list:
        res_num_str = f"{str(molnum).zfill(6)}_{str(resnum).zfill(6)}"
        summary_path = os.path.join(p_dir, f"{res_num_str}.txt")
        out_path = os.path.join(r_dir, f"{res_num_str}.out")
        
        # Check if summary exists (already completed)
        #print(f"looking for prev at {summary_path}")
        if os.path.exists(summary_path):
        #    print(f"[SKIP] {res_num_str} already completed.")
            continue
        
        # Check if output exists but summary doesn't (process it)
        if os.path.exists(out_path):
            # Check if output is valid before processing
            try:
                with open(out_path) as f:
                    text = f.read()
                    if "completed on" in text.lower():
                        print(f"[PROCESS] Valid output exists for {res_num_str}, processing...")
                        process_result(out_path, molnum, resnum, p_dir, native_lambda)
                        continue
                    else:
                        print(f"[WARNING] Output exists but incomplete for {res_num_str}, will rerun...")
                        #os.remove(out_path)  # Remove invalid output
            except Exception:
                print(f"[WARNING] Could not read output for {res_num_str}, will rerun...")
                #os.remove(out_path)
        
        # Prepare new input files
        prepared = prepare_residue_files(molnum, resnum, mae_path, in_path, r_dir)
        if prepared:
            res_num_str, in_copy_path = prepared
            in_files.append(f"{res_num_str}.in")
            print(f"Prepared {res_num_str}")
    
    if not in_files:
        print("No jobs to run!")
        return
    
    # Run jaguar on all input files
    print(f"Calculating SPE for {len(in_files)} residues")
    cmd = ['jaguar', 'run'] + in_files + ['-PARALLEL', str(int(n_cpu)), '-optimize_cpus', '-WAIT']
    print(f"Running in {r_dir}: {' '.join(cmd)}")
    
    try:
        p = subprocess.Popen(cmd, cwd=r_dir)
        p.wait()
        
        if p.returncode != 0:
            print(f"Jaguar run failed with return code {p.returncode}")
            return
        
        print("Jaguar calculations completed. Processing results...")
        
        # Process results for all residues (process_result will skip if summary exists)
        for molnum, resnum in molnum_resnum_list:
            res_num_str = f"{str(molnum).zfill(6)}_{str(resnum).zfill(6)}"
            out_path = os.path.join(r_dir, f"{res_num_str}.out")
            if process_result(out_path, molnum, resnum, p_dir, native_lambda):
                # Remove related items in r_dir if successful
                base = out_path.rsplit('.', 1)[0]
                print(f"Successfully processed {res_num_str}, removing {base}.*")
                for file_path in glob.glob(f"{base}.*"):
                    #os.remove(file_path)
                    print(f"  Removed {file_path}")
            else:
                print(f"Failed to process {res_num_str}, keeping files for debugging")

    except Exception as e:
        print(f"Error during Jaguar run: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process MAE and input files")
    parser.add_argument('-i', '--in_path', type=str, required=True, help='Path to .in file')
    parser.add_argument('-c', '--num_cpus', type=int, default=16, help='Number of CPU processors to use per process')
    parser.add_argument('-p', '--num_processes', type=int, default=4, help='Number of concurrent processes to spawn (default: 4)')
    parser.add_argument('-d','--distance',type=float, default=None,help='Limit range from QM region')
    parser.add_argument('--protein_only',action='store_true',help='Only include protein residues (excludes everything else)')
    parser.add_argument('-sasl','--screen_asl', type=str, default=None,help='Manually define residues to screen for electrostatic contributions')
    parser.add_argument('-pyp','--python_parallel',action='store_true',help='use python parallelism instead of jaguar')
    args = parser.parse_args()

    in_path = os.path.realpath(args.in_path)
    print(f"in_path: {in_path}")
    utils.textscrape.add_mae_charges_yes(in_path)
    
    with open(in_path) as f:
        in_text = f.read()
    
    native_out_path = in_path.replace(".in", ".out")
    native_lambda = utils.textscrape.extract_first_wavelength(native_out_path)
    print(f"Native lambda: {native_lambda}")

    mae_path = (re.findall(r"MAEFILE:\s*(\S+)", in_text))[0]
    spe_base = os.path.basename(in_path.rstrip('.in'))
    
    mae_st = structure.StructureReader.read(mae_path)
    qm_asl = utils.qm_connectivity.evaluate_qm_region(mae_st, in_text)
    if args.screen_asl: 
        manual_asl=args.screen_asl
    else:
        manual_asl=f'protein AND NOT {qm_asl}'
        print(manual_asl)

    molnum_resnum_list = get_nearby_mol_res(mae_st, qm_asl, args.distance, 
                                           protein_only=args.protein_only, 
                                           manual_asl=manual_asl)

    print(f"evaluating {len(molnum_resnum_list)} residues")
    
    p_dir = f"parameters_{spe_base}"
    r_dir = f"residues_{spe_base}"
    os.makedirs(p_dir, exist_ok=True)
    os.makedirs(r_dir, exist_ok=True)
    
    if args.python_parallel:
        process_all_residues(mae_path, r_dir, in_path, p_dir, args.num_processes, args.num_cpus, native_lambda, molnum_resnum_list)
    else:
        calc_jaguar_parallel(mae_path, r_dir, in_path, p_dir, args.num_processes, native_lambda, molnum_resnum_list)
        
    shutil.rmtree(r_dir)
