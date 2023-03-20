# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 17:16:38 2021

@author: android-2d
"""

import matplotlib.pyplot as plt
import numpy as np


def Existe(archivos):
    from os import path
    if path.exists(archivos):
        # print("encontrado: ",archivos)
        return True
    else:
        # print("no encontrado: ",archivos)
        return False


def GraficarBarras(inicio=1, fin=100, conservados=[]):
    for g, matrix in conservados.items():
        x = matrix[0]
        largo = len(x)
        y = 100
        x_labels = matrix[1]
        plt.figure(figsize=(14, 6))
        plt.bar(x, y)
        plt.xlabel("Aminoacidos")
        plt.ylabel("Porcentaje de conservación")
        titulo = "Análisis del virus del Papiloma Humano protenia " + str(g)
        plt.title(titulo)
        plt.xticks(x, x_labels)
        # Es muy largo para graficar todo solo muestro una parte
        plt.xlim(inicio, fin)
        salida = "IMG/" + str(g) + ".png"
        plt.savefig(salida)
        plt.show()


def RegenerarMatriz(matrix):
    a = np.array(matrix[0]).astype(int)
    largo = int(max(a)) + 1
    x = np.array(range(largo))
    y = np.array(range(largo)).astype(str)
    y[:] = ""
    for i in range(len(matrix[0])):
        posicion = int(matrix[0][i])
        label = matrix[1][i]
        x[posicion] = posicion
        y[posicion] = label
    matrix_full = np.array([x, y])
    return largo, matrix_full


def MenuRango(largo, riesgo, proteina):
    print()
    print("El alineamiento ", riesgo, "Proteina: ", proteina, "posee una longitud de: ", largo)
    print("Que desea realizar:")
    print("1- Graficar la totalidad de la longitud")
    print("2- Graficar un rango")
    print("3- Volver al menu anterior")
    salir = False
    inicio = 0
    fin = 0
    opcion = int(input("Ingrese una opción: "))
    if opcion == 1:
        inicio = 0
        fin = largo
    if opcion == 2:
        inicio = int(input("Ingrese Posición de inicio: "))
        fin = int(input("Ingrese Posición de fin: "))
    if opcion == 3:
        salir = True
    return salir, inicio, fin


def ListarGenomas():
    import CargaBD
    c = 1
    for nombre, k in CargaBD.genomas.items():
        print("Genoma número ", c, ": ", nombre)
        c += 1
    print()


def ListarProteinas():
    import CargaBD
    c = 1
    for nombre, k in CargaBD.proteinas.items():
        print("Proteina número ", c, ": ", nombre, " - Cepa: ", str(nombre[:-2]), " - ID: ", k)
        c += 1
    print()


def GraficarBarrasPares2(leyenda_proteina="", matrix=0, x=[], y=[], riesgo_cepa=""):
    x_labels = x
    x = [count for count, value in enumerate(x)]
    plt.figure(figsize=(14, 6))
    plt.bar(x, y)
    plt.xlabel("Cepas")
    plt.ylabel("Bit score")
    titulo = "Comparaciones modelos HMM de alto riesgo vs " + riesgo_cepa + "riesgo - Proteina: " + leyenda_proteina
    plt.title(titulo)
    plt.xticks(x, x_labels)
    salida = "IMG/" + riesgo_cepa + leyenda_proteina + ".png"
    plt.savefig(salida)
    plt.show()
    print("Se ha guardado el gráfico en: ", salida)
    print()


def GraficarBarrasPares(g="", matrix=0, riesgo="", proteina="", porcentaje=90):
    a = np.array(matrix[0]).astype(int)
    largo = int(max(a)) + 1
    salir, inicio, fin = MenuRango(largo, riesgo, proteina)
    leyenda_proteina = " Proteina:"
    if proteina == "":
        leyenda_proteina = " "
    if salir == False:
        x = np.array(matrix[0])
        y = np.array(matrix[1])
        for count, valor in enumerate(y):
            if valor != "":
                y[count] = porcentaje
            else:
                y[count] = 0
        y = y.astype(int)
        x_labels = matrix[1]
        plt.figure(figsize=(14, 6))
        plt.bar(x, y)
        plt.xlabel("Aminoacidos")
        plt.ylabel("Porcentaje de conservación")
        titulo = "Análisis del virus del Papiloma Humano - Clasificación: " + riesgo + leyenda_proteina + proteina
        plt.title(titulo)
        plt.xticks(x, x_labels)
        # Es muy largo para graficar todo solo muestro una parte
        plt.xlim(inicio, fin)
        salida = "IMG/" + str(g) + ".png"
        plt.savefig(salida)
        plt.show()
        print("Se ha guardado el gráfico en: ", salida)
        print()


def GraficarBarrasV1(g="", matrixFiltrada=0, riesgo="", proteina="", porcentaje=90):
    largo, matrix = RegenerarMatriz(matrixFiltrada)
    salir, inicio, fin = MenuRango(largo, riesgo, proteina)
    leyenda_proteina = " Proteina:"
    if proteina == "":
        leyenda_proteina = " "
    if salir == False:
        x = np.array(matrix[0])
        y = np.array(matrix[1])
        for count, valor in enumerate(y):
            if valor != "":
                y[count] = porcentaje
            else:
                y[count] = 0
        y = y.astype(int)
        x_labels = matrix[1]
        plt.figure(figsize=(14, 6))
        plt.bar(x, y)
        plt.xlabel("Aminoacidos")
        plt.ylabel("Porcentaje de conservación")
        titulo = "Análisis del virus del Papiloma Humano - Clasificación: " + riesgo + leyenda_proteina + proteina
        plt.title(titulo)
        plt.xticks(x, x_labels)
        # Es muy largo para graficar todo solo muestro una parte
        plt.xlim(inicio, fin)
        salida = "IMG/" + str(g) + ".png"
        plt.savefig(salida)
        plt.show()
        print("Se ha guardado el gráfico en: ", salida)
        print()


def GraficarMatriz(matriz, genoma1, genoma2):
    fig, ax = plt.subplots()
    plt.xlabel(genoma1)
    plt.ylabel(genoma2)
    titulo = "Análisis del virus del Papiloma Humano " + genoma1 + " vs " + genoma2
    plt.title(titulo, fontsize=8)
    plt.figure(figsize=(16, 10))
    im = ax.imshow(matriz, interpolation='nearest', origin='lower')
    fig.colorbar(im)
    salida = "IMG/" + genoma1 + "vs" + genoma2 + ".png"
    plt.savefig(salida)
    plt.show()
    print("Se ha guardado el gráfico en: ", salida)
    # plt.imshow(Z2, cmap ="Greens", alpha = 0.7, interpolation ='bilinear', extent = extent)
    # plt.savefig(titulo+'.png', bbox_inches='tight')
    # plt.close()


def GraficarMatriz2(matriz, genoma1, genoma2):
    fig, ax = plt.subplots()
    plt.xlabel(genoma1)
    plt.ylabel(genoma2)
    titulo = "Análisis del virus del Papiloma Humano "
    plt.title(titulo, fontsize=8)
    plt.figure(figsize=(16, 10))
    im = ax.imshow(matriz, interpolation='nearest', origin='lower')
    ax.set_xticklabels(genoma1)

    fig.colorbar(im)
    salida = "IMG/.png"
    plt.savefig(salida)
    plt.show()
    print("Se ha guardado el gráfico en: ", salida)
    # plt.imshow(Z2, cmap ="Greens", alpha = 0.7, interpolation ='bilinear', extent = extent)
    # plt.savefig(titulo+'.png', bbox_inches='tight')
    # plt.close()


"""
def ExisteG(archivos):
    from os import remove
    from os import path
    if path.exists(archivos):
        print("encontrado: ",archivos)
        return True    
"""


def Main(parametros):
    global hola
