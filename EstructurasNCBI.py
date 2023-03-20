#!/usr/bin/env python
# coding: utf-8

# # TP6 Grupo Nº 5
# Perez Ernesto Rafael rafaelperezctes@gmail.com Sofia Erdozain sofierdozain@gmail.com
# Created on Wed Nov  3 14:02:59 2021

# @author: Grupo Nº 5


import sys
import Alineamientos
import CargaBD
import Main as principal
import Utils

proteinas = CargaBD.proteinas
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
from io import StringIO

def buscarBD(archivo):
    f_record = next(SeqIO.parse(archivo, "fasta"))
    # busco en la BD
    print("corriendo blast recuperando datos...")
    result_handle = NCBIWWW.qblast("blastp", "pdb", f_record.format("fasta"))
    # guardo resultados
    qblast_out = 'm_cold_blast.out'
    archivo_p = "BD/" + qblast_out
    with open(archivo_p, "w") as save_file:
        blast_results = result_handle.read()
        save_file.write(blast_results)
    result_handle.close()
    return blast_results

def AnalizarResultados(blast_results, nombre):
    print("Analizando los resultados y extrayendo información...")
    # opción 1 - abre el archivo guardado para analizarlo
    # opción 2 - crear un identificador a partir de la cadena y analizarlo
    string_result_handle = StringIO(blast_results)
    b_record = NCBIXML.read(string_result_handle)
    # ahora obtengo la información de alineación
    # para todos los valores e mayores que algún umbral
    E_VALUE_THRESH = 0.5
    lista_blast = []
    print("Se encontraron ", len(b_record.alignments), "estructuras asociadas a ", nombre)
    c = 1
    for alignment in b_record.alignments:
        for hsp in alignment.hsps:
            # if hsp.expect < E_VALUE_THRESH:
            print()
            print("Proteina ", c, " ID:", alignment.accession)
            lista_blast.append(str(alignment.accession))
            print("****Alineacion****")
            print("secuencia: %s" % alignment.title)
            print("longitud: %i" % alignment.length)
            print("e valor: %f" % hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
        c += 1
    return lista_blast


def SubMenu2(lista_blast):
    print()
    print("### Modulo: Estructuras NCBI ###")
    print()
    indice = int(input("Elija el numero de la estructura que desea descargar - 0 para salir: "))
    if indice == 0:
        SubMenu1()
    else:
        from Bio.PDB import PDBList
        pdbl = PDBList()
        lista_pdb = [lista_blast[indice - 1][0:4]]  # debo buscar solo las 4 letras iniciales
        pdbl.download_pdb_files(lista_pdb, obsolete=False, pdir="BD/PDB/", file_format="pdb", overwrite=True)
        # pdbl.retrieve_pdb_file('1ms5', file_format="pdb")
        # https://files.rcsb.org/download/1MS5.pdb


def SubMenu1():
    print("### Modulo: Estructuras NCBI ###")
    print()
    Utils.ListarProteinas()
    proteinas_lista = list(CargaBD.proteinas)
    indice = int(input("Elija el numero de la proteina para analizar: "))
    nombre = str(proteinas_lista[indice - 1])
    archivo = "BD/" + nombre + ".fasta"
    return archivo, nombre


def Visualizar():
    import subprocess
    from subprocess import PIPE
    print(
        "Se ejecutara Jupyter-netbook en su navegador y se abrira una consola (Cierre la consola de Jupyter para volver al menu )")
    comando = "jupyter notebook EstructurasNCBI-visualizador.ipynb"
    result = subprocess.Popen(comando)
    result.communicate()
    return 1


def Menu():
    print()
    print("### Modulo: Estructuras NCBI###")
    print()
    print("Que desea realizar:")
    print("1- Buscar estructura para Proteinas almacenadas")
    print("2- Visualizar proiteinas")
    print("3- Alinear (Beta)")
    print("4- Volver al menu principal")
    print("5- Salir")
    CargaBD.LeerInicio()
    opcion_principal = int(input("Ingrese una opción: "))
    if opcion_principal == 1:
        archivo, nombre = SubMenu1()
        blast_results = buscarBD(archivo)
        lista_blast = AnalizarResultados(blast_results, nombre)
        SubMenu2(lista_blast)
        Menu()
    if opcion_principal == 2:
        Visualizar()
        Menu()
    if opcion_principal == 3:
        Utils.ListarGenomas()
        Menu()
    if opcion_principal == 4:
        principal.Menu()
    if opcion_principal == 5:
        print("Gracias por utilizar BIOG5")
        sys.exit()
# Menu()
