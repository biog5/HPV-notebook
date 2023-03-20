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
mmcif_dict = LeerCIF()
    
# PDB
resolution = structure.header["resolution"]

"""
Keys: 
name, head, deposition_date, release_date, structure_method, resolution, structure_reference (which maps to a list of references), journal_reference, author, compound (which maps to a dictionary with various information about the crystallized compound), has_missing_residues, and missing_residues.
"""

model = structure[0]
print(dir(model))
cadena=((model.get_chains))
m=(structure.get_chains())
chain = model['A']
residue = chain[100]
residue.get_resname()    # returns the residue name, e.g. "ASN"
residue.is_disordered()  # returns 1 if the residue has disordered atoms
residue.get_segid()      # returns the SEGID, e.g. "CHN1"
print(residue.id)
print(residue.full_id)
a = structure[0]["A"][100]["CA"]
print(a.full_id)
a.get_name()       # atom name (spaces stripped, e.g. "CA")
a.get_id()         # id (equals atom name)
a.get_coord()      # atomic coordinates
a.get_vector()     # atomic coordinates as Vector object
a.get_bfactor()    # isotropic B factor
a.get_occupancy()  # occupancy
a.get_altloc()     # alternative location specifier
a.get_sigatm()     # standard deviation of atomic parameters
a.get_siguij()     # standard deviation of anisotropic B factor
a.get_anisou()     # anisotropic B factor
a.get_fullname()   # atom name (with spaces, e.g. ".CA.")
print(dir(residue))
print(residue.id[0])

# grupo HEMO    
# 1 a-i
residues = model.get_residues()
for residue in residues:
    if (residue.get_resname() == "HEM"):
        print(residue.get_id())
        print(residue.get_full_id())
        for atom in residue:
            print(atom.get_fullname() )
            
# 1 a-ii           
r_no_name= mmcif_dict["_pdbx_entity_nonpoly.name"]
r_no_id= mmcif_dict["_pdbx_entity_nonpoly.comp_id"]
residues = model.get_residues()
for residue in residues:
    if (residue.get_resname() == "SO4") or (residue.get_resname() == "OXY" ):
        print(residue.get_id())
        print(residue.get_full_id())



# 1 b
structure = LeerPDB()
cadenas=0
residuos=0
atomos=0
atomos_agua=0
atomos_proteina=0
residuos_proteina =0
residuos_hetero = 0
atomos_hetero = 0 
print(atom)
for model in structure:
    for chain in model:
        # cadenas totales
        cadenas +=1
        for residue in chain:
            #residuos totales
            residuos +=1
            #atomos totales
            for atom in residue:
                atomos +=1
                if (atom.get_id() == "N") :
                   print(residue.get_full_id())
            #atomos de proteinas
            if (residue.id[0] == " "):
                residuos_proteina +=1
                for atom in residue:
                    atomos_proteina +=1
            else: #heteroatomos
                residuos_hetero +=1
                for atom in residue:
                    atomos_hetero+=1                   
            if (residue.id[0] == "W"):
                for atom in residue:
                    atomos_agua +=1
                 
print(cadenas,residuos,atomos,atomos_agua,atomos_hetero)

from Bio.PDB import Selection
atom_list = Selection.unfold_entities(chain, "A")
residue_list = Selection.unfold_entities(atom_list, "R")
chain_list = Selection.unfold_entities(atom_list, "C")


from Bio.PDB import PDBParser, CaPPBuilder
polypeptide_builder = CaPPBuilder()
counter = 1
for polypeptide in polypeptide_builder.build_peptides(structure):
    seq = polypeptide.get_sequence()
    print(f"Sequence: {counter}, Length: {len(seq)}")
    print(seq)
    counter += 1

from Bio.SeqUtils.ProtParam import ProteinAnalysis
analyzed_seq = ProteinAnalysis(str(seq))
analyzed_seq.molecular_weight()
analyzed_seq.gravy()
analyzed_seq.count_amino_acids()
analyzed_seq.get_amino_acids_percent()
analyzed_seq.secondary_structure_fraction() # helix, turn, sheet
from Bio.SeqUtils.ProtParam import ProtParamData
analyzed_seq.protein_scale(window=7, param_dict=ProtParamData.kd)

# Create empty list for chains
all_seqs = []
counter = 1
# For each polypeptide in the structure, run protein analysis methods and store in dict
for pp in polypeptide_builder.build_peptides(structure):
    seq_info = {} # create an empty dict
    seq = pp.get_sequence() # get the sequence like above
    analyzed_seq = ProteinAnalysis(str(seq)) # needs to be a str
# Specify dict keys and values    
    seq_info['Sequence Number'] = counter # set sequence id
    seq_info['Sequence'] = seq # store BioPython Seq() object
    seq_info['Sequence Length'] = len(seq) # length of seq
    seq_info['Molecular Weight'] = analyzed_seq.molecular_weight()
    seq_info['GRAVY'] = analyzed_seq.gravy() # hydrophobicity 
    seq_info['AA Count'] = analyzed_seq.count_amino_acids() 
    seq_info['AA Percent'] = analyzed_seq.get_amino_acids_percent()
# tuple of (helix, turn, sheet)
    seq_info['Secondary Structure'] = \
        analyzed_seq.secondary_structure_fraction()
    
    # Update all_seqs list and increase counter
    all_seqs.append(seq_info)
    counter += 1
all_seqs[0] # select the first sequence

"""
from Bio.PDB import *
import nglview as nv
import ipywidgets
view = nv.show_biopython(structure)
"""
from Bio.PDB import *
########
a1 = structure[0]["A"][100]["CA"]
a2 = structure[0]["A"][101]["CA"]
a3 = structure[0]["A"][102]["CA"]
a4 = structure[0]["A"][103]["CA"]

distance = a1-a2

vector1 = a1.get_vector()
vector2 = a2.get_vector()
vector3 = a3.get_vector()
angle = calc_angle(vector1, vector2, vector3)
print(distance,angle)

vector4 = a4.get_vector()
#Y si son 4 podemos medir el diedro...
angle = calc_dihedral(vector1, vector2, vector3, vector4)
print(angle)

listname=[]
residues = model.get_residues()
for residue in residues:
    #print(residue.get_resname())
    listname.append(residue.get_resname())
    if (residue.get_resname() == "H0H"):
        print(residue.get_id())
        print(residue.get_full_id())
        r=residue
        for atom in residue:
            print(atom.get_fullname() )


print(structure[0]["A"][("H_HEM",155," ")]["FE"].get_coord())
x = (structure[0]["A"][("H_HEM",155," ")]["O2A"])
#print(structure[0]["A"][("H_HEM",155," ")]["N"].get_coord())
o= (structure[0]["A"][("H_OXY",555," ")]["O1"].get_coord())
print(o)
residues = model.get_residues()


#n = residues['N'].get_vector()
# basado en doc de biopython
from Bio.PDB import PDBIO, Select
class CASelect(Select):
 def accept_atom(self, atom):
     if atom.get_name()=='CA':
         return True
     else:
         return False

io = PDBIO()
io.set_structure(structure)
io.save('CA_only.pdb', CASelect())


import nglview as nv
structure = LeerPDB()
view = nv.show_biopython(structure)
view
view2 = nv.show_biopython(structure)
view2.clear_representations()
view2.add_representation('cartoon', selection='protein', color='blue')
view2.add_representation('licorice', selection='not hydrogen')
view2



from Bio.PDB import Selection
def LeerPDB2():
    parser = PDBParser(PERMISSIVE=1)
    structure_id1 = 'miPDB'  
    filename1 = 'p38A.pdb'
    structure1 = parser.get_structure(structure_id1, filename1)
    
    structure_id2 = 'miPDB'  
    filename2 = 'p38B.pdb'
    structure2 = parser.get_structure(structure_id2, filename2)    
    return structure1, structure2

structure1, structure2 = LeerPDB2()

# Obtener id de modelo y cadena
for model in structure1:
    print(model.get_id())
    for chain in model:
        print(chain.get_id())
                
model_referencia = structure1[0]
chain_referencia = model_referencia['A']
atom_referencia = Selection.unfold_entities(chain_referencia, "A")

model_movil = structure2[0]
chain_movil = model_movil['A']
atom_movil = Selection.unfold_entities(chain_movil, "A")

#verificacion en caso de tener distinta longitud
if len(atom_referencia) > len(atom_movil):    
    diferencia= len(atom_referencia) - len(atom_movil)
    atom_referencia=atom_referencia[:-diferencia]
elif len(atom_referencia) < len(atom_movil):
    diferencia= len(atom_movil) - len(atom_referencia)
    atom_movil=atom_movil[:-diferencia]
    
#calculo la matriz
sup = Superimposer()
sup.set_atoms(atom_referencia, atom_movil)

#Puedo imprimir la matriz y el RMSD
print(sup.rotran)
print(sup.rms)

# Aplico la matriz sobre la estructura MOVIL
sup.apply(atom_movil)


import numpy

def DistanciaR1R2(residue_one, residue_two) :
    atomo1="CA" # debo elegir uno para poder buscar en ambos
    atomo2="CA"
    #atomo1="CB" 
    #atomo2="CB"
    atomo1="O" 
    atomo2="N"
    diff_vector  = residue_one[atomo1].coord - residue_two[atomo1].coord
    distancia = numpy.sqrt(numpy.sum(diff_vector * diff_vector))
    return distancia

def MatrixDistancia(chain_one, chain_two) :
    answer = numpy.zeros((len(chain_one), len(chain_two)), numpy.float)
    for fila, residue_one in enumerate(chain_one) :
        for columna, residue_two in enumerate(chain_two) :
            answer[fila, columna] = DistanciaR1R2(residue_one, residue_two)
    return answer

matriz = MatrixDistancia(chain_referencia, chain_movil)
mapa_contacto = matriz < 8 # referencia de la guia

print ("Minima distancia", numpy.min(matriz)) # cer creo el profe dijo menor a 4 sino estan unidos
print ("Maxima distancia", numpy.max(matriz))

# graficar
import numpy as np
import matplotlib.pyplot as plt
plt.imshow(matriz, cmap='gray', label="distancia en A") # probar leyenda
plt.colorbar()
plt.title("Mapa de contacto CA")
plt.ylabel("Cadena 1")
plt.xlabel("Cadena 2") 
plt.show()

# superficie
# se debe instalar el promgram msms
# conda install msms
# conda update msms

from Bio.PDB.ResidueDepth import ResidueDepth
from Bio.PDB.ResidueDepth import get_surface
def LeerPDB2():
    parser = PDBParser()
    structure_id = 'miPDB'  
    filename = '1mbo.pdb'
    structure = parser.get_structure(structure_id, filename)
    return structure
structure=LeerPDB2()

import Bio.PDB.ResidueDepth as res_depth
parser = PDBParser()
structure = parser.get_structure('miPDB',"1mbo.pdb")
model = structure[0]
rd = ResidueDepth(model, '1mbo.pdb') #1
rd = res_depth(model)  # 2
surface = rd.surface

import Bio.PDB.ResidueDepth as res_depth
from Bio.PDB.ResidueDepth import min_dist
structure = LeerPDB()
model = structure[0] # 3
surface = get_surface(model)
#distancia a la superficie 
coord = (1.113, 35.393,  9.268) #´Prueba
dist = min_dist(coord, surface)

# Distancia promedio de todos los atomos de un residuo
chain = model['A']
res152 = chain[152]
rd = res_depth(res152, surface)
rd2 = ResidueDepth(model, '1mbo.pdb')

residues = model.get_residues()
for residue in residues:  
    print(res_depth(residue, surface))
    residue_depth, ca_depth = rd2[residue]


#conda install -c salilab dssp
#conda install -c salilab -c speleo3 dssp 
#conda install -c salilab dssp
#https://biopython.org/docs/1.76/api/Bio.PDB.DSSP.html
#https://swift.cmbi.umcn.nl/gv/dssp/

"""
Bio.PDB.ResidueDepth.ca_depth( residuo , superficie )
Devuelve la profundidad de CA.

claseBio.PDB.ResidueDepth.ResidueDepth( modelo , pdb_file = Ninguno )
Bases: Bio.PDB.AbstractPropertyMap.AbstractPropertyMap

"""
from Bio.PDB import *

from Bio.PDB.DSSP import DSSP
dssp = DSSP(model,"1mbo.pdb", dssp='mkdssp')



