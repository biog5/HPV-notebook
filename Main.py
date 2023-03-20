#!/usr/bin/env python
# coding: utf-8

# # TP6 Grupo Nº 5
# Perez Ernesto Rafael rafaelperezctes@gmail.com Sofia Erdozain sofierdozain@gmail.com
# Created on Wed Nov  3 14:02:59 2021

# @author: Grupo Nº 5

import sys


def Menu():
    print("#####################################")
    print("###### Bienvenidos a BIOG5 ######")
    print("#####################################")
    print()
    print("Que desea realizar:")
    print("1- Alineamientos de a pares usando Blast")
    print("2- Alineamientos multiples usando Clustal Omega")
    print("3- Modelos HMM")
    print("4- Analisis de variantes en proteinas")
    print("5- Analisis de variantes en genomas")
    print("6- Estructuras NCBI")
    print("7- Configuar Riesgos")
    print("8- Salir")
    opcion_principal = int(input("Ingrese una opción: "))

    if opcion_principal == 1:
        import Alineamientos
        Alineamientos.Menu()

    if opcion_principal == 2:
        import AlineamientosMSA
        AlineamientosMSA.Menu()

    if opcion_principal == 3:
        import ModelosHMM
        ModelosHMM.Menu()

    if opcion_principal == 4:
        import VariantesProteinas
        VariantesProteinas.Menu()

    if opcion_principal == 5:
        import VariantesGenomas
        VariantesGenomas.Menu()

    if opcion_principal == 6:
        import EstructurasNCBI
        EstructurasNCBI.Menu()

    if opcion_principal == 7:
        import ConfigurarRiesgos
        ConfigurarRiesgos.Menu()

    if opcion_principal == 8:
        print("Gracias por utilizar BIOG5")
        sys.exit()


if __name__ == '__main__':
    Menu()
