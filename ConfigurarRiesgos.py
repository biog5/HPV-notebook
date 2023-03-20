#!/usr/bin/env python
# coding: utf-8

# # TP6 Grupo Nº 5
# Perez Ernesto Rafael rafaelperezctes@gmail.com Sofia Erdozain sofierdozain@gmail.com
# Created on Wed Nov  3 14:02:59 2021

# @author: Grupo Nº 5

from Bio import Entrez
import Utils
import sys
import CargaBD
import Main as principal
import csv

genomas = {}
unspecified = None
high_risk = None
low_risk = None
proteinas = {}
leyendas = ["Alto Riesgo", "Bajo Riesgo", "No Determinado Riesgo"]
CargaBD.LeerInicio()


def LeerListaClasificacionGuardar(bd='clasificacion.csv'):
    archivo = "BD/" + bd
    global unspecified, high_risk, low_risk
    handle = open(archivo)
    reader = list(csv.reader(handle))
    lista = {}
    for row in reader:
        # print(row)
        lista[row[0]] = row[1:]
    handle.close()
    unspecified = lista["unspecified"]
    high_risk = lista["high_risk"]
    low_risk = lista["low_risk"]


LeerListaClasificacionGuardar()


def LeerListaClasificacionA():
    bd = 'clasificacion.csv'
    archivo = "BD/" + bd
    with open(archivo, newline='') as File:
        reader = csv.reader(File)
        for count_a, row in enumerate(reader):
            print("Lita de ", row[0], " riesgo")
            for count_b, value in enumerate(row[1:]):
                print(count_b + 1, ":", value)


def ListaIndividual(row, leyenda):
    print("Lista de cepas de", leyenda)
    for count_b, value in enumerate(row):
        print(count_b + 1, ":", value)


def ListaIndividualN(archivo, columna_a, columna_b, leyenda):
    print("Lista de cepas de", leyenda)
    handle = open(archivo)
    reader = list(csv.DictReader(handle))
    # print("Leyendo lista de proteina, espere por favor...")
    c = 1
    for row in reader:
        key = row[columna_a]
        valor = row[columna_b]
        print(c, key, valor)
        c += 1
    handle.close()
    return reader


def Alta(row, valor):
    row.append(valor)
    return row


def Baja(row, indice):
    leyenda = "Confirma eliminar el elemento " + str(row[indice]) + " si - no: "
    confirmacion = input(leyenda)
    if confirmacion == "si" or confirmacion == "s":
        row.pop(indice)
    else:
        print("Tarea cancelada")
    return row


def Modificacion(row, valor, indice):
    row[indice] = valor
    return row


def RestaurarG():
    bd_g0 = 'genomas-backup.csv'
    archivo_g0 = "BD/" + bd_g0
    bd_g1 = 'genomas.csv'
    archivo_g1 = "BD/" + bd_g1
    handle_g0 = open(archivo_g0)
    reader_g0 = list(csv.DictReader(handle_g0))
    with open(archivo_g1, 'w', newline='') as csvfile:
        fieldnames = ['#Genoma', 'IdGenoma']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in reader_g0:
            writer.writerow(row)
    handle_g0.close()


def RestaurarP():
    bd_g0 = 'proteinas-backup.csv'
    archivo_g0 = "BD/" + bd_g0
    bd_g1 = 'proteinas.csv'
    archivo_g1 = "BD/" + bd_g1
    handle_g0 = open(archivo_g0)
    reader_g0 = list(csv.DictReader(handle_g0))
    with open(archivo_g1, 'w', newline='') as csvfile:
        fieldnames = ['#Genoma', 'IdProtein']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in reader_g0:
            writer.writerow(row)
    handle_g0.close()


def GrabarClasificacion():
    import csv
    global unspecified, high_risk, low_risk
    bd = 'clasificacion.csv'
    archivo = "BD/" + bd
    with open(archivo, 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        spamwriter.writerow(["unspecified"] + unspecified)
        spamwriter.writerow(["low_risk"] + low_risk)
        spamwriter.writerow(["high_risk"] + high_risk)
    print("Writing complete")


def SubMenuBaja(row):
    print()
    indice = int(input("Ingrese una indice del elemento (0 Volver atras): "))
    if indice == 0:
        Menu()
    else:
        elemento = row[indice - 1]
        salida = Baja(row, indice - 1)
        print("Se ha eliminado correctamente el item", elemento)
    return salida


def SubMenuAlta(row):
    print()
    valor = input("Ingrese la nueva cepa formato HPV'numero' (0 Volver atras): ")
    if valor == "0":
        Menu()
    else:
        salida = Alta(row, valor)
        print("Se ha agregado correctamente el item", valor)
    return salida


def AltaN(row, genoma, ide, columna_a, columna_b):
    dicc = {columna_a: genoma, columna_b: ide}
    row.append(dicc)
    return row


def SubMenuNAP(row, columna_a, columna_b):
    print()
    cepa = input("Ingrese la nueva cepa formato HPV'numero' (0 Volver atras): ")
    proteina = input("Ingrese la proteina correspondiente (E1,E2,E7,L1,L2)' (0 Volver atras): ")
    ide = input("Ingrese el identificador NCBI ejemplo (CAA75467.1, NP_040287.1, etc)' (0 Volver atras): ")
    if cepa == "0" or proteina == "0" or cepa == "0":
        Menu()
    else:
        genoma = cepa + proteina
        salida = AltaN(row, genoma, ide, columna_a, columna_b)
        print("Se ha agregado correctamente el item", genoma, ide)
    return salida


def SubMenuNAG(row, columna_a, columna_b):
    print()
    genoma = input("Ingrese la nueva cepa formato HPV'numero' (0 Volver atras): ")
    ide = input("Ingrese el identificador ref_seq NCBI ejemplo (CAA75467.1, NP_040287.1, etc)' (0 Volver atras): ")
    if genoma == "0" or ide == "0":
        Menu()
    else:
        salida = AltaN(row, genoma, ide, columna_a, columna_b)
        print("Se ha agregado correctamente el item", genoma, ide)
    return salida


def SubMenuModificacion(row):
    print("### Modulo: Configurar Riesgos ###")
    print()
    indice = int(input("Ingrese una indice del elemento (0 Volver atras): "))
    if indice == 0:
        Menu()
    else:
        valor = input("Ingrese la nueva cepa formato HPV'numero' (0 Volver atras): ")
        if valor != "0":
            valor_anterior = row[indice - 1]
            salida = Modificacion(row, valor, indice - 1)
            print("Se ha modificado correctamente", valor_anterior, "por", valor)
    return salida


def SubMenuModificacionN(row, columna_a, columna_b):
    print("### Modulo: Configurar Riesgos ###")
    print()
    indice = int(input("Ingrese una indice del elemento (0 Volver atras): "))
    if indice == 0:
        Menu()
    else:
        genoma = input("Ingrese la nueva cepa formato HPV'numero' (0 Volver atras): ")
        ide = input(
            "Ingrese Ingrese el identificador ref_seq NCBI ejemplo (CAA75467.1, NP_040287.1, etc) ' (0 Volver atras): ")
        if genoma != "0" and ide != "0":
            valor = {columna_a: genoma, columna_b: ide}
            valor_anterior = row[indice - 1]
            salida = Modificacion(row, valor, indice - 1)
            print("Se ha modificado correctamente", valor_anterior, "por", valor)
    return salida


def SubMenuModificacionNP(row, columna_a, columna_b):
    salida = []
    print("### Modulo: Configurar Riesgos ###")
    print()
    indice = int(input("Ingrese una indice del elemento (0 Volver atras): "))
    if indice == 0:
        Menu()
    else:
        cepa = input("Ingrese la nueva cepa formato HPV'numero' (0 Volver atras): ")
        proteina = input("Ingrese la proteina correspondiente (E1,E2,E7,L1,L2)' (0 Volver atras): ")
        ide = input("Ingrese el identificador NCBI ejemplo (CAA75467.1, NP_040287.1, etc)' (0 Volver atras): ")
    if cepa != "0" or proteina != "0" or cepa != "0":
        genoma = cepa + proteina
        valor = {columna_a: genoma, columna_b: ide}
        valor_anterior = row[indice - 1]
        salida = Modificacion(row, valor, indice - 1)
        print("Se ha modificado correctamente", valor_anterior, "por", valor)
    return salida


def GuardarGlobal(lista, indice):
    global unspecified, high_risk, low_risk
    if indice == 0:
        high_risk = lista
    if indice == 1:
        low_risk = lista
    if indice == 2:
        unspecified = lista
    GrabarClasificacion()


def RestaurarC():
    LeerListaClasificacionGuardar('clasificacion-backup.csv')
    GrabarClasificacion()
    print("Se han restaurado correctamente las listas de riesgos")


def GrabarDicc(lista_nueva, archivo='proteinas.csv', columna_a="#Genoma", columna_b="IdProtein", leyenda="Proteinas"):
    import csv
    with open(archivo, 'w+', newline='') as csvfile:
        fieldnames = [columna_a, columna_b]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for fila in lista_nueva:
            writer.writerow(fila)


def SubMenu1():
    global unspecified, high_risk, low_risk
    print()
    operacion_leyenda = ["Agregar", "Modificar", "Eliminar"]
    print("### Modulo: Configurar Riesgos ###")
    print("Podra editar y reclasificar las cepas conocidas y disponibles")
    print("Que desea realizar sobre los item:")
    print("1- Agregar")
    print("2- Modificar")
    print("3- Eliminar")
    print("4- Restaurar configuracion ")
    print("5- Volver al menu anterior")
    print("6- Volver al menu principal")
    print("7- Salir")
    LeerListaClasificacionGuardar()
    operacion = int(input("Ingrese una opción: "))
    if operacion == 1:
        lista, indice = SubMenu2(operacion_leyenda[operacion - 1])
        ListaIndividual(lista, leyendas[indice])
        lista_nueva = SubMenuAlta(lista)
        GuardarGlobal(lista_nueva, indice)
        ListaIndividual(lista_nueva, leyendas[indice])
    if operacion == 2:
        lista, indice = SubMenu2(operacion_leyenda[operacion - 1])
        ListaIndividual(lista, leyendas[indice])
        lista_nueva = SubMenuModificacion(lista)
        GuardarGlobal(lista_nueva, indice)
        ListaIndividual(lista_nueva, leyendas[indice])
    if operacion == 3:
        lista, indice = SubMenu2(operacion_leyenda[operacion - 1])
        ListaIndividual(lista, leyendas[indice])
        lista_nueva = SubMenuBaja(lista)
        GuardarGlobal(lista_nueva, indice)
        ListaIndividual(lista_nueva, leyendas[indice])
    if operacion == 4:
        RestaurarC()
    if operacion == 5:
        Menu()
    if operacion == 6:
        principal.Menu()
    if operacion == 7:
        print("Gracias por utilizar BIOG5")
        sys.exit()


def SubMenuN1():
    global unspecified, high_risk, low_risk
    operacion_leyenda = ["Agregar", "Modificar", "Eliminar"]
    CargaBD.LeerInicio()
    proteina_lista = CargaBD.proteinas_dictReader
    genoma_lista = CargaBD.genomas_dictReader
    bd_p = 'proteinas.csv'
    archivo_p = "BD/" + bd_p
    columna_a_p = "#Genoma"
    columna_b_p = "IdProtein"
    bd_g = 'genomas.csv'
    archivo_g = "BD/" + bd_g
    columna_a_g = "#Genoma"
    columna_b_g = "IdGenoma"
    print()
    print("### Modulo: Configurar Riesgos ###")
    print("Podra editar identificadores de genomas y proteinas")
    print("Que desea realizar sobre los item:")
    print("1- Agregar genoma")
    print("2- Modificar genoma")
    print("3- Eliminar genoma")
    print("4- Agregar proteina")
    print("5- Modificar proteina")
    print("6- Eliminar proteina")
    print("7- Restaurar configuracion genomas")
    print("8- Restaurar configuracion proteinas")
    print("9- Volver al menu anterior")
    print("10- Volver al menu principal")
    print("11- Salir")
    operacion = int(input("Ingrese una opción: "))

    if operacion == 1:
        genoma_lista = ListaIndividualN(archivo_g, columna_a_g, columna_b_g, "Genomas")
        lista_nueva = SubMenuNAG(genoma_lista, columna_a_g, columna_b_g)
        GrabarDicc(lista_nueva, archivo_g, columna_a_g, columna_b_g, "Genomas")
        ListaIndividualN(archivo_g, columna_a_g, columna_b_g, "Genomas")

    if operacion == 2:
        genoma_lista = ListaIndividualN(archivo_g, columna_a_g, columna_b_g, "Genomas")
        lista_nueva = SubMenuModificacionN(genoma_lista, columna_a_g, columna_b_g)
        GrabarDicc(lista_nueva, archivo_g, columna_a_g, columna_b_g, "Genomas")
        ListaIndividualN(archivo_g, columna_a_g, columna_b_g, "Genomas")

    if operacion == 3:
        genoma_lista = ListaIndividualN(archivo_g, columna_a_g, columna_b_g, "Genomas")
        lista_nueva = SubMenuBaja(genoma_lista)
        GrabarDicc(lista_nueva, archivo_g, columna_a_g, columna_b_g, "Genomas")
        ListaIndividualN(archivo_g, columna_a_g, columna_b_g, "Genomas")

    if operacion == 4:
        proteina_lista = ListaIndividualN(archivo_p, columna_a_p, columna_b_p, "Proteinas")
        lista_nueva = SubMenuNAP(proteina_lista, columna_a_p, columna_b_p)
        GrabarDicc(lista_nueva, archivo_p, columna_a_p, columna_b_p, "Proteinas")
        ListaIndividualN(archivo_p, columna_a_p, columna_b_p, "Proteinas")

    if operacion == 5:
        proteina_lista = ListaIndividualN(archivo_p, columna_a_p, columna_b_p, "Proteinas")
        lista_nueva = SubMenuModificacionNP(proteina_lista, columna_a_p, columna_b_p)
        GrabarDicc(lista_nueva, archivo_p, columna_a_p, columna_b_p, "Proteinas")
        ListaIndividualN(archivo_p, columna_a_p, columna_b_p, "Proteinas")

    if operacion == 6:
        proteina_lista = ListaIndividualN(archivo_p, columna_a_p, columna_b_p, "Proteinas")
        lista_nueva = SubMenuBaja(proteina_lista)
        GrabarDicc(lista_nueva, archivo_p, columna_a_p, columna_b_p, "Proteinas")
        ListaIndividualN(archivo_p, columna_a_p, columna_b_p, "Proteinas")

    if operacion == 7:
        RestaurarG()
    if operacion == 8:
        RestaurarP()
    if operacion == 9:
        Menu()
    if operacion == 10:
        principal.Menu()
    if operacion == 11:
        print("Gracias por utilizar BIOG5")
        sys.exit()


def SubMenu2(leyenda):
    print()
    print("Sobre que lista desea", leyenda, " elementos:")
    print("1- lista de alto riesgo")
    print("2- lista de bajo riesgo")
    print("3- lista de no determinado riesgo")
    LeerListaClasificacionGuardar()
    global unspecified, high_risk, low_risk, leyendas
    opcion_principal = int(input("Ingrese una opción: "))
    if opcion_principal == 1:
        lista = high_risk
        indice = 0
    if opcion_principal == 2:
        lista = low_risk
        indice = 1
    if opcion_principal == 3:
        lista = unspecified
        indice = 2
    return lista, indice


def Menu():
    print()
    print("### Modulo: Configurar Riesgos ###")
    print()
    print("Que desea realizar:")
    print("1- Listar riesgos almacenados")
    print("2- Configuar Conocidas")
    print("3- Configuar Nuevas")
    print("4- Volver al menu principal")
    print("5- Salir")
    opcion_principal = int(input("Ingrese una opción: "))
    if opcion_principal == 1:
        LeerListaClasificacionA()
        Menu()
    if opcion_principal == 2:
        SubMenu1()
        Menu()
    if opcion_principal == 3:
        SubMenuN1()
        Menu()
    if opcion_principal == 4:
        principal.Menu()
    if opcion_principal == 5:
        print("Gracias por utilizar BIOG5")
        sys.exit()
# Menu()
