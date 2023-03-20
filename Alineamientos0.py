#!/usr/bin/env python
# coding: utf-8

# # TP6 Grupo Nº 5
# Perez Ernesto Rafael rafaelperezctes@gmail.com Sofia Erdozain sofierdozain@gmail.com
# Created on Wed Nov  3 14:02:59 2021

# @author: Grupo Nº 5


import Bio.SeqIO as bsio
import matplotlib.pyplot as plt
import numpy as np
import Utils
import Main as principal
import sys
import CargaBD
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

warnings.filterwarnings("ignore")


# %matplotlib inline
#  https://www.papilocare.com/copia-de-caracteristicas

def LeerBD():
    if not (('unspecified_risk' and "low_risk" and "high_risk") in globals()):
        global unspecified_risk
        global low_risk
        global high_risk
        CargaBD.LeerInicio()
        unspecified_risk = CargaBD.unspecified
        low_risk = CargaBD.low_risk
        high_risk = CargaBD.high_risk


LeerBD()

############# 2: Alineamientos Blast
# Compararemos cada protenia (e1,e2,e7,l1,l2) de cada cepa de alto riesgo
# contra cada protenia (e1,e2,e7,l1,l2) de cepas de low and unspecified riesgo
# Para esto dividiremos en 3 partes: 
# 2.1-Agrupación  
# 2.2-blast
# 2.3-Analisis

archivos_temp = []
archivos_agrupados = CargaBD.archivos_agrupados
archivos_agrupados_gen = CargaBD.archivos_agrupados_gen
archivos_MSA = []
bandera_agrupada = 0

############################# 2.1: Agrupacion
# Leer
genes = ["E1", "E2", "E7", "L1", "L2"]
control_grabados = []


def LeerYGrabar(proteina, salida):
    global control_grabados
    grabacion = [str(proteina) + str(salida)]
    if grabacion not in control_grabados:
        control_grabados.append([str(proteina) + str(salida)])
        archivo = 'BD/' + proteina + '.fasta'
        proteinas = list(bsio.parse(archivo, 'fasta'))
        bsio.write(proteinas, salida, 'fasta')


def BuscarNombres(genoma, gen):
    genoma = genoma.upper()
    for nombre, k in CargaBD.proteinas.items():
        # tamanio1=len(nombre)
        # tamanio2=len(genoma+gen)+1
        intermedio = nombre[len(genoma):-(len(gen))]
        if (genoma in nombre) and (gen in nombre):
            if len(intermedio) <= 1 and (intermedio.isdigit() == False):
                return nombre


def Limpiar(archivos):
    from os import remove
    from os import path
    for i in archivos:
        if path.exists(i):
            # print("Eliminando: ",i)
            remove(i)


def AgruparRiesgos():
    # Primero grabo las protenias ojetivo de alto riesgo
    global archivos_agrupados, bandera_agrupada
    if bandera_agrupada == 0:
        archivos_agrupados = []  # Aunque lo tengo, lo uso por si actualizo el codigo y necesito mas clasificaciones
        for genoma in high_risk:  # nombre del genoma
            for gen in genes:
                salida_nombre = 'BD/' + "high_risk" + gen + '.fasta'
                salida_nombre_gen = 'BD/' + gen + '.fasta'
                archivo = BuscarNombres(str(genoma), str(gen))
                if archivo != None:
                    archivo_salida = open(salida_nombre, "a+")
                    archivo_salida_gen = open(salida_nombre_gen, "a+")
                    LeerYGrabar(archivo, archivo_salida)
                    LeerYGrabar(archivo, archivo_salida_gen)
                    if salida_nombre not in archivos_agrupados:
                        archivos_agrupados.append(salida_nombre)
                    if salida_nombre_gen not in archivos_agrupados_gen:
                        archivos_agrupados_gen.append(salida_nombre_gen)
                    archivo_salida.close()
                    archivo_salida_gen.close()

                    # Segundo grabo las protenias ojetivo de bajo riesgo
        for genoma in low_risk:  # nombre del genoma
            for gen in genes:
                salida_nombre = 'BD/' + "low_risk" + gen + '.fasta'
                salida_nombre_gen = 'BD/' + gen + '.fasta'
                archivo = BuscarNombres(str(genoma), str(gen))
                if archivo != None:
                    archivo_salida = open(salida_nombre, "a+")
                    archivo_salida_gen = open(salida_nombre_gen, "a+")
                    LeerYGrabar(archivo, archivo_salida)
                    LeerYGrabar(archivo, archivo_salida_gen)
                    if salida_nombre not in archivos_agrupados:
                        archivos_agrupados.append(salida_nombre)
                    if salida_nombre_gen not in archivos_agrupados_gen:
                        archivos_agrupados_gen.append(salida_nombre_gen)
                    archivo_salida.close()
                    archivo_salida_gen.close()

                    # Tercero grabo las protenias ojetivo de no especificado riesgo
        for genoma in unspecified_risk:  # nombre del genoma
            for gen in genes:
                salida_nombre = 'BD/' + "unspecified_risk" + gen + '.fasta'
                salida_nombre_gen = 'BD/' + gen + '.fasta'
                archivo = BuscarNombres(str(genoma), str(gen))
                if archivo != None:
                    archivo_salida = open(salida_nombre, "a+")
                    archivo_salida_gen = open(salida_nombre_gen, "a+")
                    LeerYGrabar(archivo, archivo_salida)
                    LeerYGrabar(archivo, archivo_salida_gen)
                    if salida_nombre not in archivos_agrupados:
                        archivos_agrupados.append(salida_nombre)
                    if salida_nombre_gen not in archivos_agrupados_gen:
                        archivos_agrupados_gen.append(salida_nombre_gen)
                    archivo_salida.close()
                    archivo_salida_gen.close()
                ############################# 2.2: Blastp
    bandera_agrupada = 1


import subprocess
from subprocess import PIPE, Popen
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
from io import StringIO
import Bio.SearchIO as bpio

"""
Intalar Blast Windows:
guia: https://www.ncbi.nlm.nih.gov/books/NBK52637/
release: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.0/
configurar la variable de entorno:BLASTDB_LMDB_MAP_SIZE=1000000

Ejemplo parametros
-perc_identity 75 #identidad
-qcov_hsp_perc 55 #cobertura querry 
-outfmt "7  #formato
"""

# globales
comparaciones_full = [{}, {}, {}]
detalles_aln = {}
dicc_querry = {}
cont_querry = 0
contador_comparaciones = 0
Limpiar(archivos_agrupados)  # Solo una vez al iniciar y lo dejo xq lo necesita MSA
Limpiar(archivos_agrupados_gen)
Limpiar(archivos_temp)
control_grabados = []


def reinciarBlast():
    global archivos_agrupados, archivos_agrupados_gen, bandera_agrupada, comparaciones_full, archivos_temp
    comparaciones_full = [{}, {}, {}]
    detalles_aln = {}
    dicc_querry = {}
    cont_querry = 0
    contador_comparaciones = 0
    archivos_temp = []
    archivos_agrupados = CargaBD.archivos_agrupados
    archivos_agrupados_gen = CargaBD.archivos_agrupados_gen
    bandera_agrupada = 0


def Recorrer(resultados, genoma, gen, riesgo):
    c = 0
    lista_aux = []
    lista_hit = []
    global dicc_querry
    global cont_querry
    dicc_hsp = {}
    cont_hsp = 0
    for query in resultados:
        cepa_query = str(query.description)[-4:-1].replace(" ", "")
        # solo corro el primero que blast tre ordenado por el mejor
        if not (cepa_query[0].isdigit()):
            cepa_query = cepa_query[1]
        cepa_query = "HPV" + str(cepa_query).upper()

        for hit in query:
            parte_numero = (str(hit.description)[-4:-1]).replace(" ", "")
            if not (parte_numero[0].isdigit()):
                parte_numero = parte_numero[1]
            cepa2 = "HPV" + str(parte_numero).upper()
            # 1 hacer uno solo
            # agregar id=int cont -- quizas un diccionario
            cont_hsp = 0
            dicc_hsp = {}
            for hsp in hit:
                # print(cepa2, hsp.bitscore)
                # if cepa2=="HPV68a":
                # opcion_principal = int(input("Ingrese una opción: "))
                # por cada hit hay varias diagonales
                # puedo elegir la mejor bitscore del hit y mostrar detalles
                # 2 guardar datos de los hsp con
                # hsp.ident_num
                # hsp.pos_num
                # hsp.gap_num
                # hsp.aln_span
                if cepa_query == "HPV7":
                    cepa_query = "HPV18"
                if cepa2 == "HPV7":
                    cepa2 = "HPV18"
                indice = str(riesgo) + cepa2 + gen + cepa_query
                dicc_hsp[cont_hsp] = [hsp.bitscore, hsp.ident_num, hsp.pos_num, hsp.gap_num, hsp.aln_span, cepa_query,
                                      cepa2, genoma, gen]
                lista_aux.append([cepa2, int(hsp.bitscore)])
                c += 1
                cont_hsp += 1
            dicc_querry[indice] = dicc_hsp
            cont_querry += 1
            # break # para elegir un solo hit el primero
    return cepa_query, lista_aux


def AlinearRiesgos():
    global cont_querry
    cont_querry = 0
    for gen in genes:
        archivo_low = "BD/low_risk" + gen + '.fasta'
        archivo_unspecified = "BD/unspecified_risk" + gen + '.fasta'
        archivo_high2 = "BD/high_risk" + gen + '.fasta'

        for genoma in high_risk:  # por cada genoma y gen de alto rieago
            archivo_entrada = BuscarNombres(str(genoma), str(gen))
            if archivo_entrada != None:
                archivo_high = 'BD/' + archivo_entrada + '.fasta'
                archivo_salida = "BD/" + genoma + gen + "_high_vs_" + gen + "_low.blast"
                archivo_salida2 = "BD/" + genoma + gen + "_high_vs_" + gen + "_unspecified.blast"
                archivo_salida3 = "BD/" + genoma + gen + "_high_vs_" + gen + "_high.blast"

                result = subprocess.run(
                    "blastp -query " + archivo_high + " -subject " + archivo_low + " -out " + archivo_salida,
                    stdout=PIPE)
                print("Corriendo BLASTP ---> Cepa:", genoma, " Proteina:", gen, " Contra grupo: Bajo riego ")
                result2 = subprocess.run(
                    "blastp -query " + archivo_high + " -subject " + archivo_unspecified + " -out " + archivo_salida2,
                    stdout=PIPE)
                print("Corriendo BLASTP ---> Cepa:", genoma, " Proteina:", gen, " Contra grupo: No especificado riego ")

                result3 = subprocess.run(
                    "blastp -query " + archivo_high + " -subject " + archivo_high2 + " -out " + archivo_salida3,
                    stdout=PIPE)
                print("Corriendo BLASTP ---> Cepa:", genoma, " Proteina:", gen, " Contra grupo: Alto riego ")

                archivos_temp.append(archivo_salida)
                archivos_temp.append(archivo_salida2)
                archivos_temp.append(archivo_salida3)
                resultados = list(bpio.parse(archivo_salida, 'blast-text'))
                resultados2 = list(bpio.parse(archivo_salida2, 'blast-text'))
                resultados3 = list(bpio.parse(archivo_salida3, 'blast-text'))

                cepa_query, lista_aux = Recorrer(resultados, genoma, gen, 0)
                cepa_query2, lista_aux2 = Recorrer(resultados2, genoma, gen, 1)
                cepa_query3, lista_aux3 = Recorrer(resultados3, genoma, gen, 2)
                comparaciones_full[0][archivo_entrada] = lista_aux
                comparaciones_full[1][archivo_entrada] = lista_aux2
                comparaciones_full[2][archivo_entrada] = lista_aux3

    print("Calculos terminados se analizaran resultados: ")


############################# 2.3: Análisis
lista_cepas = []
tamanio = 17
matriz = np.zeros((tamanio, tamanio))
import pandas as pd


def GraficarA(objetivo, x, y, proteina, cepa):
    plt.figure(figsize=(12, 6))
    print("Graficando y guardando analisis", " Cepa: ", cepa, " Proteina: ", proteina)
    plt.bar(x, y)
    plt.xlabel("Cepas de bajo riesgo")
    plt.ylabel("BitScore")
    titulo = "Análisis del virus del Papiloma Humano - " + "Cepa: " + cepa + " Proteina: " + proteina
    titulo2 = "Análisis del virus del Papiloma Humano - " + "Cepa " + cepa + " Proteina " + proteina
    plt.title(titulo)
    archivo_grafica = "IMG/" + titulo2 + ".png"
    plt.savefig(archivo_grafica)


dicc_proteina = {}
dicc_aux = {}


def AnalizarB(leyenda, E1, E2, E7, L1, L2, count=0, dicc_valores=0):
    global contador_comparaciones, dicc_aux
    global lista_cepas, matriz

    for i in ['E1', 'E2', 'E7', 'L1', 'L2']:
        print()
        print("Relaciones evolutivas con mejores scores encontradas para la proteina", i)
        print("Cepa", leyenda, "Riesgo ______ Cepa Alto Riesgo")
        x = eval(i).items()
        for versus, altos in eval(i).items():
            print(leyenda, ":", versus)
            # print("altos:",altos)
            # debor re ordenar altos
            altos_aux = []
            for j in altos:
                indice = str(count) + str(versus) + i + str(j)
                datos = dicc_querry.get(indice)
                bitscore = datos[0][0]
                altos_aux.append([j, bitscore])
            altos_aux = sorted(altos_aux, key=lambda x: x[1], reverse=True)
            # altos= [x[0] for x in altos_aux ]
            dicc_aux[versus] = altos[:4]
            for j in altos:
                # indice=dicc_valores.get(str(count)+str(versus)+i+str(j))
                indice = str(count) + str(versus) + i + str(j)
                datos = dicc_querry.get(indice)
                bitscore = datos[0][0]
                positivos = datos[0][1]
                gaps = datos[0][3]
                print("|______indice = ", contador_comparaciones, "cepa:", j, " bitscore:", bitscore, " positivos:",
                      positivos, " gaps:", gaps)
                contador_comparaciones += 1
        dicc_proteina[(i, count)] = dicc_aux
        dicc_aux = {}
    # comentario


def InicializarDiccionarios():
    E1 = {}
    E2 = {}
    E7 = {}
    L1 = {}
    L2 = {}
    return E1, E2, E7, L1, L2


def AnalizarA():
    E1, E2, E7, L1, L2 = InicializarDiccionarios()
    leyendas = ["Bajo", "No especificado", "Alto"]
    dicc_valores = {}
    global dicc_aux
    c = 0
    for count, comparaciones in enumerate(comparaciones_full):
        for objetivo, k in comparaciones.items():
            if len(k) > 0:
                matriz_aux = np.array(k)
                eje_x = matriz_aux[:, 0]
                eje_y = list(map(int, matriz_aux[:, 1]))
                proteina = str(objetivo[-2:])
                cepa = str(objetivo[:-2])
                eje_x0 = eje_x[0]
                cepa_baja = str(eje_x0)
                if cepa == cepa_baja:
                    eje_x0 = eje_x[1]
                valor = eval(proteina).get(eje_x0)
                if valor != None:
                    valor.append(cepa)
                else:
                    valor = [cepa]
                eval(proteina)[eje_x0] = valor
                # GraficarA(objetivo,eje_x,eje_y, proteina, cepa)
                dicc_valores[str(count) + cepa_baja + proteina + str(valor[-1])] = c
                c += 1
        AnalizarB(leyendas[count], E1, E2, E7, L1, L2, count, dicc_valores)
        E1, E2, E7, L1, L2 = InicializarDiccionarios()
        valor = []
        # dicc_aux = {}


# Limpiar(archivos_agrupados) # Solo una vez al iniciar y lo dejo xq lo necesita MSA
# Limpiar(archivos_temp)


def Main():
    print("1: Obteniendo datos y agrupandolos por riesgo, espere por favor ...")
    # Limpiar(archivos_agrupados)  # Solo una vez al iniciar y lo dejo xq lo necesita MSA
    # Limpiar(archivos_agrupados_gen)
    Limpiar(archivos_temp)
    # Revisar que los archivos no se esten regrabando xq darian error
    AgruparRiesgos()
    #####parametrizar aqui#######
    print("2: Se alineran genes-proteinas contra grupos de riesgo, espere por favor ...")
    AlinearRiesgos()
    # aqui mostar el score
    AnalizarA()
    # aqui agregar menu de detalles incluir un codido id=cont y mostrar datos hsp
    # Utils.GraficarMatriz2(matriz, lista_cepas, "genoma2")


def GraficarMatrizProteinas():
    global dicc_proteina
    dicc_proteina_final = {}
    for i in ['E1', 'E2', 'E7', 'L1', 'L2']:
        dicc_temp = {}
        dicc_temp.update(dicc_proteina[(i, 0)])
        dicc_temp.update(dicc_proteina[(i, 1)])
        dicc_temp.update(dicc_proteina[(i, 2)])
        dicc_proteina_final[i] = dicc_temp

    for clave_m, e in dicc_proteina_final.items():
        lista_filas = []
        lista_columnas = []
        matriz2 = []
        # obtengo el tamaño
        for clave, valor in e.items():
            if clave not in lista_filas:
                lista_filas.append(clave)
            for i in valor:
                if i not in lista_columnas:
                    lista_columnas.append(i)
        matriz2 = np.zeros((len(lista_filas), len(lista_columnas)))

        # cargo
        for clave, valor in e.items():
            indice_fila = lista_filas.index(clave)
            incremento = 1
            for i in valor:
                indice_columna = lista_columnas.index(i)
                matriz2[indice_fila][indice_columna] += incremento
                incremento -= 0.1

        df = pd.DataFrame(data=matriz2, index=lista_filas, columns=lista_columnas)
        sns.set()
        cmap = sns.diverging_palette(0, 230, 90, 60, as_cmap=True)
        s = sns.heatmap(df, annot=False, cmap='Blues', linewidth=0.3)  # , cmap='coolwarm', square=True
        plt.yticks(rotation=0, fontsize=9)
        plt.xticks(rotation=0, fontsize=9)
        titulo = "Mejores relaciones por proteina: " + clave_m.upper()
        s.set(xlabel="Cepas de alto riesgo", title=titulo, ylabel="Cepas HPV")
        plt.gcf().set_size_inches(21, 14)
        salida = "BD/" + str(clave_m.upper()) + "matriz.jpg"
        plt.savefig(salida)
        plt.show()
        print("Se ha guardado el gráfico en: ", salida)


def SubMenu1():
    print()
    print("### Analisis de resultados ###")
    print()
    print("Que desea realizar:")
    print("1- Graficar mejores relaciones")
    print("2- Volver al menu anterior")
    print("3- Volver al menu principal")
    print("4- Salir")
    CargaBD.LeerInicio()
    opcion_principal = int(input("Ingrese una opción: "))
    if opcion_principal == 1:
        GraficarMatrizProteinas()
        SubMenu1()
    if opcion_principal == 2:
        Menu()
    if opcion_principal == 3:
        principal.Menu()
    if opcion_principal == 4:
        print("Gracias por utilizar BIOG5")
        sys.exit()


def Menu():
    print()
    print("### Modulo: Alineamiento de a pares ###")
    print()
    print("Que desea realizar:")
    print("1- Listar Proteinas almacenadas")
    print("2- Listar Genomas almacenadas")
    print("3- Realizar alineamiento de a pares de proteinas con BLAST y analizar ")
    print("4- Volver al menu principal")
    print("5- Salir")

    opcion_principal = int(input("Ingrese una opción: "))
    if opcion_principal == 1:
        Utils.ListarProteinas()
        Menu()
    if opcion_principal == 2:
        Utils.ListarGenomas()
        Menu()
    if opcion_principal == 3:
        Main()
        SubMenu1()
        reinciarBlast()
    if opcion_principal == 4:
        principal.Menu()
    if opcion_principal == 5:
        print("Gracias por utilizar BIOG5")
        sys.exit()
# Main()
# Menu()
