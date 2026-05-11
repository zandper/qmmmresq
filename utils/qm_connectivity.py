from schrodinger import structure
from schrodinger.structutils import analyze
import re

def find_hcaps(in_text):
        block = re.search(
            r'&qmregion(.*?)&',
            in_text,
            re.DOTALL
        ).group(1)
        # grab all integer pairs
        pairs = re.findall(r'(\d+)\s+(\d+)', block)
        if not pairs:
            return None
        # split into separate lists
        hcapqm = [int(qm) for qm, mm in pairs]
        hcapmm = [int(mm) for qm, mm in pairs]
        return hcapqm, hcapmm

def qm_connectivity(st, hcapqm_atoms, hcapmm_atoms):
    visited = set()
    stack = list(hcapqm_atoms)
    # define forbidden edges (QM → MM cuts)
    forbidden = set(zip(hcapqm_atoms, hcapmm_atoms)) | set(zip(hcapmm_atoms, hcapqm_atoms))
    while stack:
        idx = stack.pop()
        if idx in visited:
            continue
        visited.add(idx)
        atom = st.atom[idx]
        for bond in atom.bond:
            nbr = bond.otherAtom(atom)
            # skip crossing QM/MM boundary
            if (atom.index, nbr.index) in forbidden:
                continue
            if nbr.index not in visited:
                stack.append(nbr.index)
    return sorted(visited)

def evaluate_qm_region(st, in_text):
    hcaps = find_hcaps(in_text)
    if hcaps is None:
        return None
    
    hcapqm_atoms, hcapmm_atoms = hcaps
    qm_atoms = qm_connectivity(st, hcapqm_atoms, hcapmm_atoms)
    
    # Build set of unique (chain, residue_number) pairs
    residue_info = set()
    for atom_idx in qm_atoms:
        atom = st.atom[atom_idx]
        resnum = atom.resnum
        chain = atom.chain
        residue_info.add((chain, resnum))
    # Build ASL
    asl_parts = [f'(chain.name "{chain}" AND res.num {resnum})' for chain, resnum in sorted(residue_info)]
    qm_residue_asl = ' OR '.join(asl_parts)
    return qm_residue_asl

if __name__ =="__main__":
    in_path = './000001_spe.in'
    #textscrape.add_mae_charges_yes(in_path)
    with open(in_path) as f:
        in_text = f.read()
    mae_path = (re.findall(r"MAEFILE:\s*(\S+)", in_text))[0]
    mae_st = structure.StructureReader.read(mae_path)
    qasl=evaluate_qm_region(mae_st,in_text)