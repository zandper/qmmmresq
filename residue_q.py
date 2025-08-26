from schrodinger import structure
from schrodinger.application.qsite import output
#from schrodinger.test import mmshare_data_file
import os
import re
import shutil
import subprocess
from multiprocessing import process
#print("USER:", os.environ.get("USER"))
#exit()
#path ="charg"

mae_path = "/home/zap22001/charge_diff/qsite_OCPO-pH-7.4.mae"
in_path = "/home/zap22001/charge_diff/qsite_OCPO-pH-7.4.in"
mae_fname=os.path.basename(mae_path)

#get molid of molecules in qmm region
with open(in_path, 'r') as file:
    molid_list = [int(m) for m in re.findall(r'\bmolid\s*(\d+)', file.read())]
print(molid_list)

mae_st = structure.StructureReader.read(mae_path)

e_dir="energies"
os.makedirs(e_dir, exist_ok=True)
r_dir="residues"
os.makedirs(r_dir, exist_ok=True)

#iterate through all molecules in structure
for mol in mae_st.molecule:
    #Skip molecules which are in qmm region
    if mol.number in molid_list:
        print(f"molecule #{mol.number} in qmm region, skipping molecule")
    else:
#DO THIS PART ASYNCHRONOUSLY
        #iterate through all residues a molecule
        for res in mol.residue[:8]:
            copy_st = mae_st.copy() # Make copy of original structure 
            atom_list = res._atoms # Get index of all atoms in residue
            for i in atom_list: # Set partial charge of all atoms in residue to 0 
                copy_st.atom[i].partial_charge = 0
                copy_st.atom[i].solvation_charge = 0 #ASK ABOUT THE SOLVATION CHARGE

            res_num_str = str(res.resnum).zfill(5)
            res_dir=os.path.join(r_dir,res_num_str)
            os.makedirs(res_dir,exist_ok=True)
            mae_copy_path = os.path.join(res_dir,mae_fname)
            in_copy_path = mae_copy_path.replace(".mae",".in")
            shutil.copyfile(in_path, in_copy_path)
            with structure.StructureWriter(mae_copy_path) as writer: #write modified .mae file
                writer.append(copy_st)
            cmd_dir = os.path.realpath(res_dir)
            print(res_dir)
            print(os.path.basename(in_copy_path))
            print(cmd_dir)
            p = subprocess.Popen(['qsite', '-WAIT', '-PARALLEL', '1',os.path.basename(in_copy_path)], cwd=cmd_dir)
            print(str(p))
            p.wait()
            print(mae_copy_path.replace(".mae",".out"))
            result = output.QSiteOutput(mae_copy_path.replace(".mae",".out"))
            print(result.energy)
            with open(os.path.join(e_dir,f"{res_num_str}.txt"),"w") as f:
                f.write(f"{res_num_str}\t{result.energy}")
            shutil.rmtree(res_dir) 
            exit()
#DO THIS PART ASYNCHRONOUSLY

            #t = structure._AtomCollection(mae_st, atom_list)
            """
            print(t[0])
            exit()
            st_temp = mae_st.copy()
            st_temp.property = st_temp.property
            print(st_temp.property)
            atom_list= res._atoms
            print(mae_st.atom[atom_list])
            """
            #nres = build.neutralize_structure(res)
            #print(dir(atoms))
            #print(res.getCode())
            #print(res.resnum) #resnum is the mol
            #atom_list= res._atoms
            #set atom charges to 0
            #print(atom_list)
            #mae_st()
            #set charge of residue equal to 0
            
            #save .mae file
            #save .in file (needs to have matching name)
            #run calc

            #print(dir(res))

            #temp_mae_st = mae_st.neutralize residue

            

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