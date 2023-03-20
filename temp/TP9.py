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
end_id   = 4000
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
#Copyright 2007 Peter Cock, all rights reserved.
# Licenced under the GPL v2 (or later at your choice)
#
# Please see this website for details:
# http://www.warwick.ac.uk/go/peter_cock/python/protein_superposition/
import Bio.PDB
import numpy

pdb_code = "1JOY"
pdb_filename = "%s.pdb" % pdb_code
pdb_out_filename = "%s_aligned.pdb" % pdb_code
seq_str = 'MAAGVKQLADDRTLLMAGVSHDLRTPLTRIRLATEMMSEQDGYLAESINKDIEECNAIIEQFIDYLR'
use_str = '-----------RTLLMAGVSHDLRTPLTRIRLATEMMSEQDGYLAESINKDI---------------'

#Make a list of booleans which specify which residues we will consider for the RMS
#Basically I am excluding a handfull of residues at each end which are very free
#(not part of an alpha helix of a beta sheet).
use = [letter != "-" for letter in use_str]
#assert len(use) == len(seq_str)

#print "Loading PDB file %s" % pdb_filename
#structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename)

print ("Everything aligned to first model...")
ref_model = structure[0]
for alt_model in structure :
    #Build paired lists of c-alpha atoms:
    ref_atoms = []
    alt_atoms = []
    for (ref_chain, alt_chain) in zip(ref_model, alt_model) :
        for ref_res, alt_res, amino, allow in zip(ref_chain, alt_chain, seq_str, use) :
            assert ref_res.resname == alt_res.resname
            assert ref_res.id      == alt_res.id
            assert amino == Bio.PDB.Polypeptide.three_to_one(ref_res.resname)
            if allow :
                #CA = alpha carbon
                ref_atoms.append(ref_res['CA'])
                alt_atoms.append(alt_res['CA'])

    #Align these paired atom lists:
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, alt_atoms)

    if ref_model.id == alt_model.id :
        #Check for self/self get zero RMS, zero translation
        #and identity matrix for the rotation.
        assert numpy.abs(super_imposer.rms) < 0.0000001
        assert numpy.max(numpy.abs(super_imposer.rotran[1])) < 0.000001
        assert numpy.max(numpy.abs(super_imposer.rotran[0]) - numpy.identity(3)) < 0.000001
    else :
        #Update the structure by moving all the atoms in
        #this model (not just the ones used for the alignment)
        super_imposer.apply(alt_model.get_atoms())

    print ("RMS(first model, model %i) = %0.2f" % (alt_model.id, super_imposer.rms))


print ("Saving aligned structure as PDB file %s" % pdb_out_filename)
io=Bio.PDB.PDBIO()
io.set_structure(structure)
io.save(pdb_out_filename)

print ("Done")



