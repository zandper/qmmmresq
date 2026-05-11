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
    """
    Get nearby molecule/residue information around QM region.
    
    Args:
        mae_st: Schrödinger Structure object
        qm_asl: ASL string selecting QM atoms
        cutoff_angstrom: Distance cutoff in angstroms (used if manual_asl not provided)
        manual_asl: ASL string that directly specifies residues to screen (no distance filtering)
        protein_only: If True, exclude solvent and membrane, only keep protein
    
    Returns:
        List of (molecule_number, resnum) tuples for nearby residues
    """
    qm_atoms = set(analyze.evaluate_asl(mae_st, qm_asl))
    if not qm_atoms:
        raise ValueError(f"No atoms selected by ASL: '{qm_asl}'")
    
    # Select atoms either by manual ASL or by distance cutoff
    if manual_asl:
        # manual_asl directly specifies what to screen - no distance filtering
        selected_atoms = set(analyze.evaluate_asl(mae_st, manual_asl))
        if not selected_atoms:
            raise ValueError(f"No atoms selected by manual ASL: '{manual_asl}'")
        nearby_atoms = selected_atoms
    elif cutoff_angstrom is not None:
        nearby_atoms = set(measure.get_atoms_close_to_subset(mae_st, qm_atoms, cutoff_angstrom))
        nearby_atoms -= qm_atoms
    else:
        raise ValueError("Either cutoff_angstrom or manual_asl must be provided")
    
    # Filter out solvent, membrane and salt if protein_only
    if protein_only:
        solvent_atoms = set(analyze.evaluate_asl(mae_st, "solvent"))
        membrane_atoms = set(analyze.evaluate_asl(mae_st, "membrane"))
        salt_atoms=set(analyze.evaluate_asl(mae_st,'(ions) AND ((chain.name " "))'))
        nearby_atoms -= solvent_atoms
        nearby_atoms -= membrane_atoms
        nearby_atoms -= salt_atoms

    # Return unique (molecule_number, resnum) pairs
    return sorted({
        (mae_st.atom[i].molecule_number, mae_st.atom[i].resnum)
        for i in nearby_atoms
    })

def calc_single_point_residue(molnum, resnum, mae_path, r_dir, mae_fname, in_path, p_dir, n_cpu, native_lambda):

    # ---- Build file identidfiers----
    res_num_str = f"{str(molnum).zfill(6)}_{str(resnum).zfill(6)}"
    res_dir = os.path.join(r_dir, res_num_str)
    mae_copy_path = os.path.join(res_dir, mae_fname)
    out_path = mae_copy_path.replace(".mae", ".out")

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

    with structure.StructureWriter(mae_copy_path) as writer:
        writer.append(mae_st)

    # ---- Run QSite ----
    cmd = ['qsite', '-WAIT', '-HOST', 'localhost','-PARALLEL', str(int(n_cpu)),os.path.basename(in_copy_path)]
    print(cmd)
    p = subprocess.Popen(cmd, cwd=res_dir)
    p.wait()

    # ---- Post-processing ----
    try:
        result = output.QSiteOutput(out_path)
        wavelength_nm = utils.textscrape.extract_first_wavelength(out_path)
        res_contrib = wavelength_nm - native_lambda
        print(f"Total_E: {result.energy}, Residue Contribution:{res_contrib}, Wavelength: {wavelength_nm}")
        summary_path = os.path.join(p_dir, f"{res_num_str}.txt")
        with open(summary_path, "w") as f:
            f.write(f"{molnum}\t{resnum}\t{result.energy}\t{res_contrib}\t{wavelength_nm}\n")
        #remove res_dir
        shutil.rmtree(res_dir)
    except Exception as e:
        print(f"Error reading {out_path}: {e}")


def process_all_residues(mae_path, r_dir, mae_fname, in_path, p_dir, num_processes, n_cpu, native_lambda, molnum_resnum_list):
    
    args = [(molnum, resnum, mae_path, r_dir, mae_fname, in_path, p_dir, n_cpu, native_lambda) for molnum, resnum in molnum_resnum_list]
    
    print(f"Selected residues: {len(args)}")
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.starmap(calc_single_point_residue, args)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Process MAE and input files")
    parser.add_argument('-i', '--in_path', type=str, required=True, help='Path to .in file')
    parser.add_argument('-c', '--num_cpus', type=int, default=16, help='Number of CPU processors to use per process')
    parser.add_argument('-p', '--num_processes', type=int, default=4, help='Number of concurrent processes to spawn (default: 4)')
    parser.add_argument('-d','--distance',type=float, default=None,help='Limit range from QM region')
    parser.add_argument('--protein_only',action='store_true',help='Only include protein residues (excludes everything else)')
    parser.add_argument('-sasl','--screen_asl', type=str, default=None,help='Manually define residues to screen for electrostatic contributions')
    args = parser.parse_args()

    in_path = os.path.realpath(args.in_path)
    utils.textscrape.add_mae_charges_yes(in_path)
    with open(in_path) as f:
        in_text = f.read()
    native_out_path = in_path.replace((".in"),(".out"))
    native_lambda = utils.textscrape.extract_first_wavelength(native_out_path)
    print(native_lambda)

    mae_path = (re.findall(r"MAEFILE:\s*(\S+)", in_text))[0]
    mae_fname = os.path.basename(mae_path)
    spe_base = os.path.basename(in_path.rstrip('.in'))
    num_processes = args.num_processes
    n_cpu = args.num_cpus
    mae_st = structure.StructureReader.read(mae_path)
    qm_asl = utils.qm_connectivity.evaluate_qm_region(mae_st,in_text)
    molnum_resnum_list = get_nearby_mol_res(mae_st, qm_asl, args.distance, protein_only=args.protein_only, manual_asl=args.screen_asl)

    print(f"evaluating {len(molnum_resnum_list)} residues")
    

    p_dir=f"parameters_{spe_base}"
    os.makedirs(p_dir, exist_ok=True)
    r_dir=f"residues_{spe_base}"
    os.makedirs(r_dir, exist_ok=True)
    
    process_all_residues(mae_path, r_dir, mae_fname, in_path, p_dir, num_processes, n_cpu, native_lambda, molnum_resnum_list)

    shutil.rmtree(r_dir)




            

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