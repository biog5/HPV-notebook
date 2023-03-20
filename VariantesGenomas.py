#!/usr/bin/env python
# coding: utf-8

# # TP6 Grupo Nº 5
# Perez Ernesto Rafael rafaelperezctes@gmail.com Sofia Erdozain sofierdozain@gmail.com
# Created on Wed Nov  3 14:02:59 2021

# @author: Grupo Nº 5

import Bio.SeqIO as bsio
import numpy as np
import pandas as pd
from Bio import Align
import Utils
import CargaBD
import Main as principal
import sys

# Paso 1 obtener los genomas
CargaBD.LeerInicio()
unspecified_risk = CargaBD.unspecified
low_risk = CargaBD.low_risk
high_risk = CargaBD.high_risk
genomas = CargaBD.genomas


########################## Algoritmo 1: Alinear 2 genomas


def Lee2Genomas(nombre1, nombre2):
    archivo1 = 'BD/' + nombre1 + '.fasta'
    archivo2 = 'BD/' + nombre2 + '.fasta'
    if (Utils.Existe(archivo1) and Utils.Existe(archivo2)):
        align1 = bsio.read(archivo1, 'fasta')
        align2 = bsio.read(archivo2, 'fasta')
        return align1, align2


def AlinearGenomasAPares(nombre1, nombre2):
    seq1, seq2 = Lee2Genomas(nombre1, nombre2)
    aligner = Align.PairwiseAligner()
    tamanio = int(len(seq1.seq))
    matriz = np.zeros((tamanio, tamanio))
    alignments = aligner.align(seq1.seq, seq2.seq)
    # Lo paso a tablas
    x = alignments[0]
    query_alineado = list(str(x).splitlines()[0])
    target_alineado = list(str(x).splitlines()[2])
    c = 0
    for i in range(tamanio):
        if (query_alineado[i] == target_alineado[i]):
            matriz[i][i] = int(1)
            c += 1
    print("Longitud del alineamiento: ", c)
    print("Score del alineamiento: ", alignments.score)
    align_pares = pd.DataFrame([query_alineado, target_alineado])
    return matriz, align_pares


def Algoritmo1():
    global genoma1
    global genoma2
    global matriz, align_pares, align_pd
    genoma1 = "HPV16"  # ejemplo
    genoma2 = "HPV18"
    print("### Se alinearan 2 genomas ###")
    genoma1 = int(input("Ingrese el numero del genoma 1: ")) - 1
    genoma2 = int(input("Ingrese el numero del genoma 2: ")) - 1
    claves = list(CargaBD.genomas.keys())
    genoma1 = claves[genoma1]
    genoma2 = claves[genoma2]
    print()
    print("Se alineran los genomas:", genoma1, "vs genoma: ", genoma2)
    # aqui aling_pd error en el varios viene del la lectura de un archivo aqui corre
    global archivos_agrupados, genomas, archivos_genomas
    print("3: Corriendo Clustal Omega, espere por favor ...")
    # CorrerClustalOmega()
    matriz, align_pd = AlinearGenomasAPares(genoma1, genoma2)
    # matriz(posicion_conservada, aminoacido)
    matriz_ab = align_pd.to_numpy()
    matriz_a = np.array([np.arange(len(matriz_ab[0])), matriz_ab[0]])
    matriz_b = np.array([np.arange(len(matriz_ab[0])), matriz_ab[1]])
    # matriz(aminoacido, aminoacido) index=posicion
    SubMenuGraficasPares(matriz, matriz_a, matriz_b)

def AgrupacionIngresados():
    Utils.ListarGenomas()
    salida =""
    lista=[]
    claves = list(CargaBD.genomas.keys())
    while salida != "salir":
       salida=input("Ingrese el numero del genoma (para terminar escriba: salir) ")
       if salida != "salir":
             salida= int(salida)-1
             nombre = claves[salida]
             archivo = 'BD/' + nombre + '.fasta'
             lista.append(archivo)
    return lista
# Algoritmo1()

########################## Algoritmo 2: Multiples genomas

from Bio import AlignIO
import subprocess
from subprocess import PIPE
import CargaBD


def ListarArchivos():
    aux_high = []
    aux_low = []
    aux_unspecified = []
    for nombre, k in genomas.items():
        archivo = 'BD/' + str(nombre) + '.fasta'
        for h in high_risk:
            if (h in nombre):
                aux_high.append(archivo)
        for l in low_risk:
            if (l in nombre):
                aux_low.append(archivo)
        for u in unspecified_risk:
            if (u in nombre):
                aux_unspecified.append(archivo)
    archivos_genomas["high"] = aux_high
    archivos_genomas["low"] = aux_low
    archivos_genomas["unspecified"] = aux_unspecified

def AgruparArchivosIngresados(lista):
    global archivos_agrupados
    for archivo_entrada in lista:
        if Utils.Existe(archivo_entrada):
            genoma = bsio.read(archivo_entrada, 'fasta')
            archivo_salida = 'BD/' + str("genomas_selecionados") + '.fasta'
            output_file = open(archivo_salida, "a")
            output_file.write(genoma.format('fasta'))
            output_file.close()
    archivos_agrupados.append(archivo_salida)
    print()

def AgruparArchivos():
    global archivos_agrupados
    for riesgo, k in archivos_genomas.items():
        for archivo_entrada in k:
            if Utils.Existe(archivo_entrada):
                genoma = bsio.read(archivo_entrada, 'fasta')
                archivo_salida = 'BD/' + str(riesgo) + '.fasta'
                output_file = open(archivo_salida, "a")
                output_file.write(genoma.format('fasta'))
                output_file.close()
        archivos_agrupados.append(archivo_salida)
    print()


# esto es pesadisimo si uso genomas   
def CorrerClustalOmega():
    for archivo in archivos_agrupados:
        archivo_entrada = archivo
        archivo_salida = archivo.replace(".fasta", "_MSA.phylip")

        # Si se tiene instalado descomentar y correr:
        # clustalomega_cline = ClustalOmegaCommandline(infile = archivo_entrada, outfile = archivo_salida, outfmt = 'phylip', verbose = True, auto = False)
        # print(clustalomega_cline)

        # Por practicidad lo corremos portable en windows
        clustal = '"clustal-omega-1.2.2-win64/clustalo.exe"' + ' -i ' + archivo_entrada + ' -o ' + archivo_salida + ' --outfmt phylip -v --force'
        print("Corriendo clustal-omega-1.2.2-win64. sobre", archivo_entrada, "  . Espere por favor....")
        result = subprocess.run(clustal, stdout=PIPE)
        archivos_MSA.append(archivo_salida)


def LeerMSA():
    """
    pruebo con un resultado ya guardado MSA pero para todos,
    recordar descomentar en Algoritmo 2 la linea de Clustal y dejar correr
    y comenentar las dos lineas siguientes
    """
    #global archivos_MSA
    #archivos_MSA = ['BD/high_riskE1_MSA.phylip']

    for archivo in archivos_MSA:
        align = list(AlignIO.read(archivo, "phylip"))
        # for alignment in align:
        # print(alignment)
        align_pd = pd.DataFrame(align)
        # filas = align_pd.shape[0]
        msa.append(align_pd)
        return align_pd


# Paso 3 calculo de porcentajes
def CargaPorcentajes():
    porcentajes_aux = {}
    c = 0
    for align_pd in msa:
        for columna in align_pd:
            cant_hist = align_pd[columna].value_counts()
            valores = [list(cant_hist.index), list((cant_hist.values * 100) / cant_hist.sum())]
            porcentajes_aux[columna] = valores
            patron[columna] = list("x")
        porcentajes[c] = porcentajes_aux
        c += 1


def CargarConservados(porcentaje=90):
    x = []
    y = []
    x_labels = []
    for p, porcentajesP in porcentajes.items():
        for i, valor in porcentajesP.items():
            if valor[1][0] > porcentaje:
                x.append(i)
                y.append(valor[1][0])
                x_labels.append(valor[0][0])
                patron[i] = valor[0]
        matrix = np.array([x, x_labels])
        x = []
        y = []
        x_labels = []
        conservados[p] = matrix


def Imprimir():
    riesgo, opcion = ListarMSA()
    matrix = conservados.get(opcion)
    if matrix is None:
        matrix = [[0, 0], [0, 0]]
    for j in range(len(matrix[0])):
        print("Riesgo: ", riesgo, "Posición: ", matrix[0][j], "Aminoacido", matrix[1][j])
    print()
#

def GraficarBarras1(porcentaje=90):
    riesgo, opcion = ListarMSA()
    matrix = conservados.get(opcion)
    if matrix is None:
        matrix = [[0, 0], [0, 0]]
    Utils.GraficarBarrasV1("", matrix, riesgo, "")


def ListarMSA():
    print("Lista de agrupacion de genomas: ")
    c = 1
    aux_clasificacion = []
    for i in archivos_agrupados:
        clasificacion = str(i[3:-6])
        aux_clasificacion.append(clasificacion)
        print("Número:", c, "Clasificación:", clasificacion)
        c += 1
    opcion = int(input("Elija un número: ")) - 1
    riesgo = aux_clasificacion[opcion]

    """
    Número: 1 Clasificación: high
    Número: 2 Clasificación: low
    Número: 3 Clasificación: unspecified
    """
    return riesgo, opcion


def SubMenuGraficasPares(matriz, matriz_a, matriz_b):
    print("Que desea realizar:")
    print("1- Graficar ambos genomas en una matriz")
    print("2- Graficar genoma 1")
    print("3- Graficar genoma 2")
    print("4- Volver al menu principal")
    print("5- Salir")
    CargaBD.LeerInicio()
    opcion_principal = int(input("Ingrese una opción: "))
    if opcion_principal == 1:
        Utils.GraficarMatriz(matriz, genoma1, genoma2)
        #SubMenuPares()
        SubMenuGraficasPares(matriz, matriz_a, matriz_b)
    if opcion_principal == 2:
        Utils.GraficarBarrasPares("", matriz_a, "", "", 50)
        #SubMenuPares()
        SubMenuGraficasPares(matriz, matriz_a, matriz_b)
    if opcion_principal == 3:
        Utils.GraficarBarrasPares("", matriz_b, "", "", 50)
        Menu()
    if opcion_principal == 4:
        principal.Menu()
    if opcion_principal == 5:
        print("Gracias por utilizar BIOG5")
        sys.exit()

def SubMenuSeleecionGenomas():
    print("Que desea realizar:")
    print("1- Selecionar multiples genomas a alinear")
    print("2- Correr alinemiento contra todos los genomas almacenados")
    print("3- Volver al menu anterior")
    print("4- Volver al menu principal")
    print("5- Salir")
    opcion_principal = int(input("Ingrese una opción: "))
    if opcion_principal == 1:
        Utils.ListarGenomas()
        lista = AgrupacionIngresados()
        AgruparArchivosIngresados(lista)
    if opcion_principal == 2:
        ListarArchivos()
        print("Se alineran la siguientes genomas, espere por favor ...")
        AgruparArchivos()
    if opcion_principal == 3:
        Menu()
    if opcion_principal == 4:
        principal.Menu()
    if opcion_principal == 5:
        print("Gracias por utilizar BIOG5")
        sys.exit()

def SubMenu1(porcentaje=90):
    print("### Modulo: Alineamiento multiples ###")
    print()
    print("Que desea realizar:")
    print("1- Imprimir lista de aminoacidos conservados por riesgo")
    print("2- Graficar riesgo")
    print("3- Volver al menu anterior")
    print("4- Volver al menu principal")
    print("5- Salir")
    CargaBD.LeerInicio()
    opcion_principal = int(input("Ingrese una opción: "))
    if opcion_principal == 1:
        Imprimir()
        SubMenu1()
    if opcion_principal == 2:
        GraficarBarras1()
        SubMenu1()
    if opcion_principal == 3:
        Menu()
    if opcion_principal == 4:
        principal.Menu()
    if opcion_principal == 5:
        print("Gracias por utilizar BIOG5")
        sys.exit()


# Menu()

def Algoritmo2():
    print("### Se alinearan multiples genomas ### ")
    porcentaje = int(input("Ingrese el porcentaje minimo de conservación de una posición: "))
    global matriz, align_pares, align_pd, archivos_MSA, archivos_genomas, genomas, conservados
    global archivos_agrupados, msa, porcentajes, patron
    archivos_MSA = []
    archivos_genomas = {}
    genomas = CargaBD.genomas
    conservados = {}
    archivos_agrupados = []
    print("1: Obteniendo archivos de Base de Datos, espere por favor ...")
    print("2: Agrupando archivos, espere por favor ...")
    print()
    SubMenuSeleecionGenomas()
    print("3: Corriendo Clustal Omega, espere por favor ...")
    CorrerClustalOmega()
    msa = []
    align_pd = LeerMSA()
    porcentajes = {}
    patron = {}
    print("4: Cargando porcentajes, espere por favor ...")
    CargaPorcentajes()
    print("5: Filtrando porcentajes con >90% de conservación, espere por favor ...")
    CargarConservados(porcentaje=90)
    # GraficarBarras()
    print("Calculos terminados analizar resultados: ")
    print()
    SubMenu1()


# Algoritmo2()

def Menu():
    print()
    print("### Modulo: Variantes Genomas ###")
    print()
    print("Que desea realizar:")
    print("1- Listar Proteinas almacenadas")
    print("2- Listar Genomas almacenadas")
    print("3- Alinear 2 genomas")
    print("4- Alinear genomas disponibles (Beta)")
    print("5- Volver al menu principal")
    print("6- Salir")
    CargaBD.LeerInicio()
    opcion_principal = int(input("Ingrese una opción: "))
    if opcion_principal == 1:
        Utils.ListarProteinas()
        Menu()
    if opcion_principal == 2:
        Utils.ListarGenomas()
        Menu()
    if opcion_principal == 3:
        Utils.ListarGenomas()
        Algoritmo1()
        Menu()
    if opcion_principal == 4:
        #Utils.ListarGenomas()
        Algoritmo2()
        Menu()
    if opcion_principal == 5:
        principal.Menu()
    if opcion_principal == 6:
        print("Gracias por utilizar BIOG5")
        sys.exit()
# Menu()
