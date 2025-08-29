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
