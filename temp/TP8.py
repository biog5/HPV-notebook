# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 14:17:44 2021

@author: android-2d
"""

# CIF
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


def LeerPDB():
    parser = PDBParser(PERMISSIVE=1)
    structure_id = 'miPDB'  
    filename = '1mbo.pdb'
    structure = parser.get_structure(structure_id, filename)
    return structure

def LeerCIF():
    mmcif_dict = MMCIF2Dict("1mbo.cif")
    for value,key in mmcif_dict.items():
        print(value,key)
    return mmcif_dict

structure = LeerPDB()
#mmcif_dict = LeerCIF()
    
# PDB
resolution = structure.header["resolution"]
keywords = structure.header["keywords"]

"""
Keys: 
name, head, deposition_date, release_date, structure_method, resolution, structure_reference (which maps to a list of references), journal_reference, author, compound (which maps to a dictionary with various information about the crystallized compound), has_missing_residues, and missing_residues.

for model in structure:
    for chain in model:
        print(chain)
        for residue in chain:
            print(residue)
      #

print(dir(residue))
model = structure[0]
print(resolution)
cadena=((model.get_chains))
m=(structure.get_chains())
"""
import Bio.PDB
# https://gist.github.com/andersx/6354971

# Seleccione qué números de residuos desea alinear
# y ponerlos en una lista
start_id = 1
end_id   = 400
atoms_to_be_aligned = range(start_id, end_id + 1)

# Iniciar el analizador
pdb_parser = Bio.PDB.PDBParser(QUIET = True)

# Obtener las estructuras
ref_structure = pdb_parser.get_structure("reference", "4xr8-16-E6.pdb")
sample_structure = pdb_parser.get_structure("samle", "6sjv-18-E6.pdb")

# Use el primer modelo en los archivos pdb para la alineación
# Cambia el número 0 si quieres alinearte con otra estructura
ref_model    = ref_structure[0]
sample_model = sample_structure[0]

# Haz una lista de los átomos (en las estructuras) que deseas alinear.
# En este caso usamos átomos CA cuyo índice está en el rango especificado
ref_atoms = []
sample_atoms = []

# Iterar todas las cadenas en el modelo para encontrar todos los residuos
for ref_chain in ref_model:
  # Iterate of all residues in each model in order to find proper atoms
  for ref_res in ref_chain:
    # Comprobar si el número de residuo ( .get_id() ) está en la lista
    print(ref_res.get_id()[1],atoms_to_be_aligned, ref_res.id[1],ref_res)

    #if ref_res.get_id()[1] in atoms_to_be_aligned:
      # Agregar átomo de CA a la lista
      #print(ref_res.get_id()[1])
    if 'CA' in ref_res:
        ref_atoms.append(ref_res['CA'])
    #for atom in ref_res:
        #print(atom)

# Do the same for the sample structure
for sample_chain in sample_model:
  for sample_res in sample_chain:
    #if sample_res.get_id()[1] in atoms_to_be_aligned:
    if 'CA' in sample_res:
        sample_atoms.append(sample_res['CA'])

# Now we initiate the superimposer:
super_imposer = Bio.PDB.Superimposer()
super_imposer.set_atoms(ref_atoms, sample_atoms)
super_imposer.apply(sample_model.get_atoms())

# Print RMSD:
print(super_imposer.rms)

# Save the aligned version of 1UBQ.pdb
io = Bio.PDB.PDBIO()
io.set_structure(sample_structure)
io.save("16-18-e.pdb")





