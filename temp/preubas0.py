#!/usr/bin/env python
# coding: utf-8

# # TP6 Grupo Nº 5
# Perez Ernesto Rafael rafaelperezctes@gmail com Sofia Erdozain sofierdozain@gmail com
# Created on Wed Nov  3 14:02:59 2021

# @author: Grupo Nº 5


import Bio.SeqIO as bsio
import matplotlib.pyplot as plt
import numpy as np
import Utils
import Main as principal
import sys
import CargaBD
#  https://www papilocare com/copia-de-caracteristicas

import matplotlib.pyplot as plt
import Alineamientos
import CargaBD
import Main as principal
import Utils
import sys
archivos_MSA = []
archivos_agrupados = []

def ObtenerAgrupaciones():
    global archivos_agrupados
    Alineamientos.Limpiar(CargaBD.archivos_agrupados)
    Alineamientos.AgruparRiesgos()  # genera fasta que son limpiados antes
    archivos_agrupados = Alineamientos.archivos_agrupados

ObtenerAgrupaciones()
dicc_id_cepa={}
def RelacionarIdCepa():
    global dicc_id_cepa
    for archivo in archivos_agrupados:
        proteinas = list(bsio.parse(archivo, 'fasta'))
        for i in proteinas:
            len_id=len(i.name)
            id=i.name
            descripcion=i.description[len_id:]
            cepa_query = str(descripcion)[-4:-1].replace(" ", "")
            if not (cepa_query[0].isdigit()):
                cepa_query = cepa_query[1]
            cepa_query = "HPV" + str(cepa_query).upper()
            if cepa_query == "HPV7":
                cepa_query = "HPV18"
            print(id, descripcion,cepa_query)
            dicc_id_cepa[id]=cepa_query
#RelacionarIdCepa()
print()
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

matriz=[[0, 2, 1, 2, 2, 1, 3, 2, 2, 3, 2, 2, 2, 2, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 0, 0, 3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 0, 0, 0, 0, 0, 0, 0, 3, 5, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 3, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 0, 0, 0, 5, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 2, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 1, 1, 0, 2, 2, 0, 1, 0, 1, 0, 0, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 5, 2, 5, 5, 4, 1, 1, 5, 5, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 2, 1, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 0, 1, 0, 0, 0, 1, 2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0],
 [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0]]

lista_cepas=['HPV42', 'HPV16', 'HPV18', 'HPV31', 'HPV33', 'HPV35H', 'HPV39', 'HPV45', 'HPV52', 'HPV58', 'HPV59', 'HPV68A', 'HPV68B', 'HPV73', 'HPV82', 'HPV44', 'HPV56', 'HPV43', 'HPV11', 'HPV6B', 'HPV74', 'HPV67', 'HPV69', 'HPV34', 'HPV66', 'HPV53', 'HPV30', 'HPV3', 'HPV62']

df = pd.DataFrame(data=matriz, index=lista_cepas, columns=lista_cepas)
plt.figure(figsize=(7, 7))
#sns.heatmap(df, annot=True, cmap='coolwarm')
#plt.show()
k=[['HPV43', 608, 3], ['HPV11', 864, 3], ['HPV44', 669, 3], ['HPV6B', 671, 3], ['HPV74', 674, 3], ['HPV42', 744, 3]]
print(k)
k = sorted(k, key=lambda x: x[1])
k.sort(key=lambda x: x[1],reverse=True)
print(k)

e1= {'HPV43': ['HPV16', 'HPV31', 'HPV52'], 'HPV44': ['HPV39', 'HPV68A', 'HPV68B'], 'HPV11': ['HPV82']}
e1={'HPV43': ['HPV58', 'HPV16', 'HPV31'], 'HPV42': ['HPV73'], 'HPV44': ['HPV33', 'HPV58', 'HPV31'], 'HPV11': ['HPV82'], 'HPV6B': ['HPV35H', 'HPV31', 'HPV58'], 'HPV74': ['HPV52', 'HPV33', 'HPV16'], 'HPV53': ['HPV56'], 'HPV57C': ['HPV33', 'HPV35H'], 'HPV34': ['HPV68B', 'HPV39', 'HPV59'], 'HPV3': ['HPV56', 'HPV52'], 'HPV66': ['HPV31'], 'HPV62': ['HPV82'], 'HPV57B': ['HPV16', 'HPV31', 'HPV18'], 'HPV69': ['HPV58', 'HPV52', 'HPV16'], 'HPV67': ['HPV59', 'HPV39'], 'HPV40': ['HPV45'], 'HPV82': ['HPV16', 'HPV52', 'HPV31'], 'HPV73': ['HPV18', 'HPV45', 'HPV56'], 'HPV39': ['HPV82', 'HPV73'], 'HPV56': ['HPV16', 'HPV45', 'HPV33'], 'HPV59': ['HPV56'], 'HPV45': ['HPV73'], 'HPV58': ['HPV59', 'HPV39', 'HPV68A'], 'HPV68B': ['HPV73', 'HPV31', 'HPV52'], 'HPV68A': ['HPV31', 'HPV35H', 'HPV16'], 'HPV35H': ['HPV39', 'HPV56', 'HPV59']}
e2={1:e1,2:e1}


for clave_m, e1 in e2.items():
    lista_filas=[]
    lista_columnas=[]
    # obtengo el tamaño
    for clave, valor in e1.items():
        if clave not in lista_filas:
            lista_filas.append(clave)
        for i in valor:
            if i not in lista_columnas:
                lista_columnas.append(i)
    matriz2 = np.zeros((len(lista_filas), len(lista_columnas)))

    # cargo
    for clave, valor in e1.items():
        indice_fila = lista_filas.index(clave)
        incremento=1
        for i in valor:
            indice_columna = lista_columnas.index(i)
            matriz2[indice_fila][indice_columna] += incremento
            incremento -= 0.1

    df = pd.DataFrame(data=matriz2, index=lista_filas, columns=lista_columnas)
    plt.figure(figsize=(34, 14))
    #sns.heatmap(df, annot=True, cmap='coolwarm', square=True)
    sns.heatmap(df, cmap='coolwarm')
    sns.set()
    plt.xlabel("Cepas de alto riesgo")
    plt.ylabel("Cepas de Bajo/no_determindo riesgo")
    plt.title("Mejores relaciones por proteina")
    plt.show()
x=None
if x:
    print("es none")
else:
    print("no es none")