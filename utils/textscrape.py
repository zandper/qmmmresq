import re

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
