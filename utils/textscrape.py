import re

def get_molids(in_path):
    #get molid of molecules in qmm region
    with open(in_path, 'r') as file:
        molid_list = [int(m) for m in re.findall(r'\bmolid\s*(\d+)', file.read())]
    return molid_list

def extract_first_wavelength(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
    # Search for wavelength in nm
    match = re.search(r'Excitation energy =.*?([\d.]+)\s*nm', content)
    
    if match:
        wavelength_nm = float(match.group(1))
        return wavelength_nm
    else:
        return "WAVELENGTH NOT FOUND"

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
