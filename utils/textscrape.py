import re

def get_molids(in_path):
    #get molid of molecules in qmm region
    with open(in_path, 'r') as file:
        molid_list = [int(m) for m in re.findall(r'\bmolid\s*(\d+)', file.read())]
    return molid_list

def extract_first_wavelength(file_path, from_orb=86, to_orb=87):
    with open(file_path) as f:
        content = f.read()

    blocks = re.split(r'(?=Restricted Singlet Excited State\s+\d+)', content)
    exc_re = re.compile(
        r'^\s*(\d+)\s*=>\s*(\d+)\s+([-+]?\d*\.\d+(?:[eE][-+]?\d+)?)\s*$',
        re.M,
    )
    wl_re  = re.compile(r'Excitation energy\s*=.*?([\d.]+)\s*nm')
    osc_re = re.compile(r'Oscillator strength,\s*f\s*=\s*([-+]?\d*\.\d+(?:[eE][-+]?\d+)?)')

    best = (0.0, None, None)          # (rel, wavelength, osc_strength)
    for b in blocks:
        x   = re.search(r'X coeff\.(.*?)(?=Y coeff\.|\n\s*\n|\Z)', b, re.S)
        wl  = wl_re.search(b)
        osc = osc_re.search(b)
        if not x or not wl:
            continue
        coeffs = {(int(m[0]), int(m[1])): float(m[2])
                  for m in exc_re.findall(x.group(1))}
        norm = sum(abs(c) for c in coeffs.values())
        if not norm or (from_orb, to_orb) not in coeffs:
            continue
        rel = abs(coeffs[(from_orb, to_orb)]) / norm
        if rel > best[0]:
            best = (rel, float(wl.group(1)),
                    float(osc.group(1)) if osc else None)

    return best[1], best[2]           # (wavelength_nm, oscillator_strength)

def add_mae_charges_yes(in_file_path):
    # Read the file
    with open(in_file_path, "r") as f:
        lines = f.readlines()
    # Check if 'use_mae_charges=YES' already exists
    if not any("use_mae_charges=YES" in line for line in lines):
        new_lines = []
        for line in lines:
            new_lines.append(line)
            if line.strip() == "&mmkey":
                # Insert after &mmkey
                new_lines.append("use_mae_charges=YES\n")
        # Write back
        with open(in_file_path, "w") as f:
            f.writelines(new_lines)
        print('.in FILE MODIFIED, ADDED "use_mae_charges=YES" ')
