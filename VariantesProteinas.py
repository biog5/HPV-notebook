#!/usr/bin/env python
# coding: utf-8

# # TP6 Grupo Nº 5
# Perez Ernesto Rafael rafaelperezctes@gmail.com Sofia Erdozain sofierdozain@gmail.com
# Created on Wed Nov  3 14:02:59 2021

# @author: Grupo Nº 5


import subprocess
import sys
from subprocess import PIPE

import Bio.SeqIO as bsio
import numpy as np
import pandas as pd
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline

import Alineamientos
import CargaBD
import Main as principal
import Utils

# Paso 0 alinear
# paso 1 alineo genomas

# archivos_agrupados = Alineamientos.archivos_agrupados
archivos_MSA = []
archivos_proteinas = {}
clasificacionTotal = {'high_riskE1': [], 'low_riskE1': [], 'unspecified_riskE1': [], 'high_riskE2': [],
                      'low_riskE2': [], 'unspecified_riskE2': [], 'high_riskE7': [
    ], 'low_riskE7': [], 'unspecified_riskE7': [], 'high_riskL1': [], 'low_riskL1': [], 'unspecified_riskL1': [],
                      'high_riskL2': [], 'low_riskL2': [], 'unspecified_riskL2': []}

unspecified_risk = CargaBD.unspecified
low_risk = CargaBD.low_risk 
high_risk = CargaBD.high_risk
proteinas = CargaBD.proteinas


def QueClaseEs(substring):
    clases = []
    for nombre, k in proteinas.items():
        genoma = nombre[:-2]
        clases.append(genoma)
    for i in clases:
        if (substring in i):
            return i


def QueRiesgoEs():
    for i in high_risk:
        i = QueClaseEs(i)
        clasificacionRiesgos[i] = "high_risk"
    for j in low_risk:
        j = QueClaseEs(j)
        clasificacionRiesgos[j] = "low_risk"
    for k in unspecified_risk:
        k = QueClaseEs(k)
        clasificacionRiesgos[k] = "unspecified_risk"


def ListarArchivos():
    # agrupo por riesgo
    for nombre, k in proteinas.items():
        gen = nombre[-2:]
        genoma = nombre[:-2]
        archivo = 'BD/' + str(nombre) + '.fasta'
        # acumular
        QueRiesgoEs()
        # print(nombre,genoma)
        if (clasificacionRiesgos.get(genoma) != None):
            valor = clasificacionRiesgos[genoma]
            key = valor + gen
            valor = clasificacionTotal.get(key)
            if valor != None:
                valor.append(archivo)
                clasificacionTotal[key] = valor


def UnirArchivos():
    print("1: Agrupando datos por tipo de Riesgo y Proteina, espere por favor ...")
    archivos_agrupados = {}
    for clasificacion, k in clasificacionTotal.items():
        archivo_salida = 'BD/' + str(clasificacion) + '.fasta'
        if Utils.Existe(archivo_salida):
            output_file = open(archivo_salida, "w+")
        else:
            output_file = open(archivo_salida, "a")
        for archivo_entrada in k:
            if Utils.Existe(archivo_entrada):
                # print("uniendo", clasificacion, "salida", archivo_salida, archivo_entrada)
                genoma = bsio.read(archivo_entrada, 'fasta')
                output_file.write(genoma.format('fasta'))
        archivos_agrupados[clasificacion] = archivo_salida
    output_file.close()
    return archivos_agrupados


# esto es pesadisimo si uso genomas


def CorrerClustalOmega():
    archivos_MSA = []
    print("2: Corriendo alineamiento multiple con Clustal Omega, espere por favor ... ")
    for key, archivo in archivos_agrupados.items():
        archivo_entrada = archivo
        archivo_salida = 'BD/' + key + "_MSA.phylip"

        # Si se tiene instalado descomentar y correr:
        clustalomega_cline = ClustalOmegaCommandline(
            infile=archivo_entrada, outfile=archivo_salida, outfmt='phylip', verbose=True, auto=False)
        # print("Corriendo: ", clustalomega_cline)

        # Por practicidad lo corremos portable en windows
        clustal = '"clustal-omega-1.2.2-win64/clustalo.exe"' + ' -i ' + \
                  archivo_entrada + ' -o ' + archivo_salida + ' --outfmt phylip -v --force'
        result = subprocess.run(clustal, stdout=PIPE)
        archivos_MSA.append(archivo_salida)
    return archivos_MSA


def LeerMSA():
    msa = []
    alineamientos_panda = {}
    print("3: Cargando alineamientos, espere por favor ...")
    for archivo in archivos_MSA:
        align = list(AlignIO.read(archivo, "phylip"))
        align_pd = pd.DataFrame(align)
        clave = archivo[3:-7]
        # print("Cargando pandas, espere por favor ...", clave)
        alineamientos_panda[clave] = align_pd
    return alineamientos_panda


# Paso 3 calculo de porcentajes


def CargaPorcentajes():
    porcentajes_aux = {}
    alineamientos_porcentajes = {}
    print("4: Cargando porcentajes, espere por favor ...")
    for key, align_pd in alineamientos_panda.items():
        for columna in align_pd:
            # obtengo todos los elementos de una columna
            # lista_elementos = align_pd[columna].value_counts().index
            # print(lista_elementos)
            cant_hist = align_pd[columna].value_counts()
            valores = [list(cant_hist.index), list(
                (cant_hist.values * 100) / cant_hist.sum())]
            porcentajes_aux[columna] = valores
            patron[columna] = list("x")
        alineamientos_porcentajes[key] = porcentajes_aux
        porcentajes_aux = {}
    return alineamientos_porcentajes


def CargarConservados(porcentaje=90):
    x = []
    y = []
    x_labels = []
    conservados = {}
    print("5: Filtrando porcentajes con >90% de conservación, espere por favor ...")
    for p, porcentajesP in alineamientos_porcentajes.items():
        for i, valor in porcentajesP.items():
            if valor[1][0] > porcentaje:
                x.append(i)
                y.append(valor[1][0])
                x_labels.append(valor[0][0])
                # print(p, valor, i)
                patron[i] = valor[0]
        matrix = np.array([x, x_labels])
        x = []
        y = []
        x_labels = []
        conservados[p] = matrix
    return conservados


def ListarProteinas():
    print("Lista de proteinas agrupadas por riesgo con > 90% de conservación: ")
    for p, valor in CargaBD.conservados.items():
        print("Número:", p, "[", valor[0], "Proteina:", valor[1], "]")
    opcion = int(input("Elija un número: "))
    clave_conservados = CargaBD.conservados[opcion][2]
    riesgo = CargaBD.conservados[opcion][0]
    proteina = CargaBD.conservados[opcion][1]
    matrix = conservados[clave_conservados]
    return riesgo, proteina, matrix, clave_conservados


def Imprimir():
    riesgo, proteina, matrix, clave_conservados = ListarProteinas()
    for j in range(len(matrix[0])):
        print(riesgo, " Proteina:", proteina, "Posición:", matrix[0][j], "Aminoacido:", matrix[1][j])
    # Imprimir()


def GraficarBarras(porcentaje):
    riesgo, proteina, matrix, clave_conservados = ListarProteinas()
    Utils.GraficarBarrasV1(clave_conservados, matrix, riesgo, proteina, porcentaje)


def SubMenu1(porcentaje=90):
    print("### Modulo: Alineamiento multiples ###")
    print()
    print("Que desea realizar:")
    print("1- Imprimir lista de aminoacidos conservados")
    print("2- Graficar Riesgos y Aminoacidos")
    print("3- Volver al menu anterior")
    print("4- Volver al menu principal")
    print("5- Salir")
    CargaBD.LeerInicio()
    opcion_principal = int(input("Ingrese una opción: "))
    if opcion_principal == 1:
        Imprimir()
        SubMenu1()
    if opcion_principal == 2:
        GraficarBarras(porcentaje)
        SubMenu1()
    if opcion_principal == 3:
        Menu()
    if opcion_principal == 4:
        principal.Menu()

    if opcion_principal == 5:
        print("Gracias por utilizar BIOG5")
        sys.exit()


# Menu()

def Main():
    print("Se procede a buscar variantes en proteinas")
    porcentaje = int(input("Ingrese el porcentaje minimo de conservación de una posición: "))
    global clasificacionRiesgos
    clasificacionRiesgos = {}
    ListarArchivos()
    global archivos_agrupados
    Alineamientos.AgruparRiesgos()
    archivos_agrupados = UnirArchivos()
    global archivos_MSA
    archivos_MSA = CorrerClustalOmega()
    global alineamientos_panda
    alineamientos_panda = LeerMSA()
    global porcentajes
    global patron
    porcentajes = {}
    patron = {}
    global alineamientos_porcentajes
    global conservados
    alineamientos_porcentajes = CargaPorcentajes()
    conservados = CargarConservados(porcentaje)
    print("Calculos terminados analizar resultados: ")
    print()
    SubMenu1(porcentaje)


# Main()

def Menu():
    print()
    print("### Modulo: Variantes Proteinas ###")
    print()
    print("Que desea realizar:")
    print("1- Listar Proteinas almacenadas")
    print("2- Listar Genomas almacenadas")
    print("3- Realizar calculos para buscar variantes")
    print("4- Volver al menu principal")
    print("5- Salir")
    CargaBD.LeerInicio()
    opcion_principal = int(input("Ingrese una opción: "))
    if opcion_principal == 1:
        Utils.ListarProteinas()
        Menu()
    if opcion_principal == 2:
        Utils.ListarGenomas()
        Menu()
    if opcion_principal == 3:
        Main()
    if opcion_principal == 4:
        principal.Menu()
    if opcion_principal == 5:
        print("Gracias por utilizar BIOG5")
        sys.exit()
# Menu()
