#!/usr/bin/env python
# coding: utf-8

# # TP6 Grupo Nº 5
# Perez Ernesto Rafael rafaelperezctes@gmail.com Sofia Erdozain sofierdozain@gmail.com
# Created on Wed Nov  3 14:02:59 2021

# @author: Grupo Nº 5

import matplotlib.pyplot as plt
import Alineamientos
import CargaBD
import Main as principal
import Utils
import sys
import Bio.SeqIO as bsio

archivos_MSA = []
archivos_agrupados = []
archivos_agrupados_gen = []

def ObtenerAgrupaciones():
    global archivos_agrupados
    Alineamientos.Limpiar(CargaBD.archivos_agrupados)
    Alineamientos.Limpiar(CargaBD.archivos_agrupados_gen)
    #global control_grabados
    Alineamientos.bandera_agrupada = 0
    Alineamientos.control_grabados=[]
    Alineamientos.AgruparRiesgos()  # genera fasta que son limpiados antes
    archivos_agrupados = Alineamientos.archivos_agrupados
    archivos_agrupados = Alineamientos.archivos_agrupados_gen

    


############################# 3: Alineamientos MSA
"""
Alineo grupos, usando ClustalO Windows
alto riesgo E1
alto riesgo E2
alto riesgo E7
alto riesgo L1
alto riesgo L2
bajo riesgo E1
bajo riesgo E2
bajo riesgo E7
bajo riesgo L1
bajo riesgo L2
no_determinado riesgo E1
no_determinado riesgo E2
no_determinado riesgo E7
no_determinado riesgo L1
no_determinado riesgo L2

Obtengo_arboles por cada uno usando lubreria phytree

"""
# http://www.clustal.org/omega/

import subprocess
from subprocess import PIPE
dicc_id_cepa={}

def RelacionarIdCepa():
    global dicc_id_cepa
    for archivo in archivos_agrupados:
        proteinas = list(bsio.parse(archivo, 'fasta'))
        for i in proteinas:
            len_id=len(i.name)
            id=i.name
            descripcion = i.description[len_id:]
            cepa_query = str(descripcion)[-4:-1].replace(" ", "")
            if not (cepa_query[0].isdigit()):
                cepa_query = cepa_query[1]
            cepa_query = "HPV" + str(cepa_query).upper()
            if cepa_query == "HPV7":
                cepa_query = "HPV18"
            dicc_id_cepa[id]=cepa_query

def CorrerClustalOmega():
    global archivos_agrupados
    global archivos_agrupados_gen
    archivos_agrupados = archivos_agrupados + archivos_agrupados_gen
    for archivo in archivos_agrupados:
        archivo_entrada = archivo
        archivo_salida = archivo.replace(".fasta", "_MSA.phylip")
        proteina = str(archivo_entrada[-8:-6])
        riesgo = str(archivo_entrada[3:-8]).replace("_", " ")

        # Si se tiene instalado descomentar y correr:
        # clustalomega_cline = ClustalOmegaCommandline(infile = archivo_entrada, outfile = archivo_salida, outfmt = 'phylip', verbose = True, auto = False)
        print("Corriendo Clustal Omega ---> Grupo:", riesgo, " Proteina:", proteina, "espere por favor ...")
        # Por practicidad lo corremos portable en windows
        # clustal='"clustal-omega-1.2.2-win64/clustalo.exe"' + ' -i '+ archivo_entrada + ' -o '+ archivo_salida +' --infmt a2m --outfmt phylip -v --force'
        clustal = '"clustal-omega-1.2.2-win64/clustalo.exe"' + ' -i ' + archivo_entrada + ' -o ' + archivo_salida + ' --outfmt phylip -v --force'
        result = subprocess.run(clustal, stdout=PIPE)
        archivos_MSA.append(archivo_salida)
    print()


def LeerMSA():
    for archivo in archivos_MSA:
        alignments = list(AlignIO.parse(archivo, "phylip"))
        # for alignment in alignments:
        # print(alignment)
        # print()


# https://biopython.org/wiki/Phylo
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

# proteinas={"HPV16E1":"NP_041327.2", "HPV16E2":"NP_041328.1", "HPV16E7":"NP_041326.1", "HPV18E1":"NP_040312.1", "HPV18E2":"NP_040313.1", "HPV18E7":"NP_040311.1", "HPV45E1":"CAA52575.1", "HPV45E2":"CAA52576.1", "HPV45E7":"CAA52574.1", "HPV31E1":"AAA46952", "HPV31E2":"AAA46953.1", "HPV31E7":"AAA46951.1", "HPV33E1":"AAA46960.1", "HPV33E2":"AAA46961.1", "HPV33E7":"AAA46959.1", "HPV35HE1":"CAA52563.1", "HPV35HE2":"CAA52564.1", "HPV35HE7":"CAA52562.1", "HPV39E1":"AAA47052.1", "HPV39E2":"AAA47053.1", "HPV39E7":"AAA47051.1", "HPV52E1":"CAA52587.1", "HPV52E2":"CAA52588.1", "HPV52E7":"CAA52586.1", "HPV56E2":"CAA52598.1", "HPV56E7":"CAA52597.1", "HPV58E1":"BAA31847.1", "HPV58E2":"BAA31848.1", "HPV58E7":"BAA31846.1", "HPV59E1":"CAA54851", "HPV59E2":"CAA54852.1", "HPV59E7":"CAA54850.1", "HPV68AE1":"AAZ39493.1", "HPV68AE2":"AAZ39494.1", "HPV68AE7":"AAZ39492.1", "HPV68BE1":"CBY85069.1", "HPV68BE2":"CBY85070.1", "HPV68BE7":"CBY85068.1", "HPV73E1":"CAA63884.1", "HPV73E2":"CAA63885.1", "HPV73E7":"CAA63883.1", "HPV82E1":"BAA90737.1", "HPV82E2":"BAA90738.1", "HPV82E7":"BAA90736.1", "HPV23E1":"AAA79410.1", "HPV23E2":"AAA79411.1", "HPV23E7":"AAA79409.1", "HPV53E1A":"NP_597792.1","HPV53E1B":"NP_597793.1", "HPV53E2":"NP_041846.1", "HPV53E7":"NP_041845.1", "HPV66E1":"AAA79501.1", "HPV66E2":"AAA79502.1", "HPV66E7":"AAA79500.1", "HPV6BE1":"NP_040298.1", "HPV6BE2":"NP_040299.1", "HPV6BE7":"NP_040297.1","HPV11E1":"AAA46929.1", "HPV11E2":"AAA46930.1", "HPV11E7":"AAA46928.1", "HPV40E1":"CAA52569.1", "HPV40E2":"CAA52570.1", "HPV40E7":"CAA52568.1","HPV43E1":"CAF05785.1", "HPV43E2":"CAF05786.1", "HPV43E7":"CAF05784.1", "HPV44E1":"AAA79459.1", "HPV44E2":"AAA79460.1", "HPV44E7":"AAA79458.1", "HP54VE1":"NP_043290.1", "HPV54E2":"NP_043291.1", "HPV54E7":"NP_043289.1", "HPV61E1":"NP_043446.1", "HPV61E2":"NP_043447.1", "HPV61E7":"NP_043445.1", "HPV70E1":"AAC54852.1", "HPV70E2":"AAC54853.1", "HPV70E7":"AAC54851.1", "HPV72E1":"AHZ58190.1", "HPV72E2":"AHZ58191.1", "HPV72E7":"AHZ58194.1", "HPV81E1":"CAF05694.1", "HPV81E2":"CAF05695.1", "HPV81E7":"CAF05693.1", "HPV34E1":"NP_041809.1", "HPV34E2":"NP_041810.1", "HPV34E7":"NP_041808.1", "HPV57CE1":"BAF80482.1", "HPV57CE2":"BAF80483.1", "HPV57CE7":"BAF80481.1", "HPV83E1":"AAD38970.1", "HPV83E2":"AAD38971.1", "HPV83E7":"AAD38969.1", "HPV16L1":"NP_041332.2", "HPV16L2":"NP_041331.2", "HPV18L1":"NP_040317.1", "HPV18L2":"NP_040316.1", "HPV45L1":"CAA52578.1", "HPV45L2":"CAA52577.1", "HPV31L1":"AAA46956.1", "HPV31L2":"AAA46955.1", "HPV33L1":"AAA46964.1", "HPV33L2":"AAA46963.1", "HPV35HL1":"CAA52566.1", "HPV35HL2":"CAA52565.1", "HPV39L1":"AAA47056.1", "HPV39L2":"AAA47055.1", "HPV52L1":"CAA52590.1", "HPV52L2":"CAA52589.1", "HPV56L1":"CAA52600.1", "HPV56L2":"CAA52599.1", "HPV58L1":"BAA31851.1", "HPV58L2":"BAA31850.1", "HPV59L1":"CAA54856.1", "HPV59L2":"CAA54855.1", "HPV68AL1":"AAZ39498.1", "HPV68AL2":"AAZ39497.1", "HPV68BL1":"CBY85074.1", "HPV68BL2":"CBY85073.1", "HPV73L1":"CAA63887.1", "HPV73L2":"CAA63886.1", "HPV82L1":"BAA90742.1", "HPV82L2":"BAA90741.1", "HPV23L1":"AAA79414.1", "HPV23L2":"AAA79413.1", "HPV53L1":"NP_041848.1", "HPV53L2":"NP_041847.1", "HPV66L1":"AAA79505.1", "HPV66L2":"AAA79504.1", "HPV6BL1":"NP_040304.1", "HPV6BL2":"NP_040303.1", "HPV11L1":"AAA46935.1", "HPV11L2":"AAA46934.1", "HPV40L1":"CAA52572.1", "HPV40L2":"CAA52571.1","HPV42E1":"AAA47043.1", "HPV42E7":"AAA47042.1" ,"HPV42L1":"AAA47048.1", "HPV42L2":"AAA47047.1","HPV43L1":"CAF05790.1", "HPV43L2":"CAF05788.1", "HPV54L1":"NP_043294.1", "HPV54L2":"NP_043293.1", "HPV61L1":"NP_043450.1", "HPV61L2":"NP_043449.1", "HPV70L1":"AAC54857.1", "HPV70L2":"AAC54856.1", "HPV72L1":"AHZ58195.1", "HPV72L2":"AHZ58196.1", "HPV81L1":"CAF05698.1", "HPV81L2":"CAF05697.1","HPV34L1":"NP_041812.1", "HPV34L2":"NP_041811.1" ,"HPV57BL1":"AAC56600.1", "HPV57CL1":"BAF80486.1", "HPV57CL2":"BAF80485.1", "HPV83L1":"AAD38974.1", "HPV83L2":"AAD38973.1", "HPV2E1":"QLM04844.1", "HPV2E2":"QLM04841.1", "HPV2E7":"QLM04843.1", "HPV2L1":"QLM04846.1", "HPV2L2":"QLM04845.1", "HPV3E1":"CAA52471.1", "HPV3E2":"CAA52472.1", "HPV3L1":"CAA52474.1", "HPV3L2":"CAA52473.1","HPV5E1":"NP_041367.1", "HPV5E2":"NP_041368.1", "HPV5E7":"NP_041366.1", "HPV5L1":"NP_041372.1", "HPV5L2":"NP_041371.1", "HPV4E1":"NP_040891.1", "HPV4E2":"NP_040892.1", "HPV4E7":"NP_040890.1", "HPV4L1":"NP_040895.1", "HPV4L2":"NP_040894.1", "HPV7E1":"NP_041856.1", "HPV7E2":"NP_041857.1", "HPV7E7":"NP_041855.1", "HPV7L1":"NP_041859.1", "HPV7L2":"NP_041858.1", "HPV9E1":"NP_041863.1", "HPV9E2":"NP_041864.1", "HPV9E7":"NP_041862.1", "HPV9L1":"NP_041866.1", "HPV9L2":"NP_041865.1", "HPV10E1":"NP_041743.1", "HPV10E2":"NP_041744.1", "HPV10E7":"NP_041742.1", "HPV10L1":"NP_041746.1", "HPV10L2":"NP_041745.1", "HPV120E1":"CAA52498.1", "HPV12E2":"CAA52499.1", "HPV12E7":"CAA52497.1", "HPV12L1":"CAA52501.1", "HPV12L2":"CAA52500.1", "HPV13E1":"CAA44649.1", "HPV13E2":"CAA44650.1", "HPV13E7":"CAA44648.1", "HPV13L1":"CAA44654.1", "HPV13L2":"CAA44653.1", "HPV14BE1":"CAA52502.1", "HPV14BE2":"CAA52503.1", "HPV14BL1":"CAA52505.1", "HPV14BL2":"CAA52504.1", "HPV15E1":"CAA52508.1", "HPV15E2":"CAA52509.1", "HPV15E7":"CAA52507.1", "HPV15L1":"CAA52511.1", "HPV15L2":"CAA52510.1", "HPV17E1":"CAA52514.1", "HPV17E2":"CAA52515.1", "HPV17E7":"CAA52513.1", "HPV17L1":"CAA52517.1", "HPV17L2":"CAA52516.1", "HPV19E1":"CAA52520.1", "HPV19E2":"CAA52521.1", "HPV19E7":"CAA52519.1", "HPV19L1":"CAA52523.1", "HPV19L2":"CAA52522.1", "HPV20E1":"AAA79389.1", "HPV20E2":"AAA79390.1", "HPV20E7":"AAA79388.1", "HPV20L1":"AAA79393.1", "HPV20L2":"AAA79392.1", "HPV21E1":"AAA79396.1", "HPV21E2":"AAA79397.1", "HPV21E7":"AAA79395.1", "HPV21L1":"AAA79400.1", "HPV21L2":"AAA79399.1", "HPV22E1":"AAA79403.1", "HPV22E2":"AAA79404.1", "HPV22E7":"AAA79402.1", "HPV22L1":"AAA79407.1", "HPV22L2":"AAA79406.1", "HPV23E1":"AAA79410.1", "HPV23E2":"AAA79411.1", "HPV23E7":"AAA79409.1", "HPV23L1":"AAA79414.1", "HPV23L2":"AAA79413.1", "HPV24E1":"AAA79417.1", "HPV24E2":"AAA79418.1", "HPV24E7":"AAA79416.1", "HPV24L1":"AAA79421.1", "HPV24L2":"AAA79420.1", "HPV25E1":"CAA52526.1", "HPV25E2":"CAA52527.1", "HPV25E7":"CAA52525.1", "HPV25L1":"CAA52529.1", "HPV25L2":"CAA52528.1", "HPV27BE1":"BAE16265.1", "HPV27BE2":"BAE16266.1", "HPV27BE7":"BAE16264.1", "HPV27BL1":"BAE16269.1", "HPV27BL2":"BAE16268.1", "HPV28E1":"AAA79424.1", "HPV28E2":"AAA79425.1", "HPV28E7":"AAA79423.1", "HPV28L1":"AAA79428.1", "HPV28L2":"AAA79427.1", "HPV29E1":"AAA79431.1", "HPV29E2":"AAA79432.1", "HPV29E7":"AAA79430.1", "HPV29L1":"AAA79435.1", "HPV29L2":"AAA79434.1", "HPV30E1":"YP_009508156.1", "HPV30E2":"YP_009508157.1", "HPV30E7":"YP_009508155.1", "HPV30L1":"YP_009508159.1", "HPV30L2":"YP_009508158.1", "HPV32E1":"NP_041803.1", "HPV32E2":"NP_041804.1", "HPV32E7":"NP_041802.1", "HPV32L1":"NP_041806.1", "HPV32L2":"NP_041805.1", "HPV36E1":"AAA79438.1", "HPV36E2":"AAA79439.1", "HPV36E7":"AAA79437.1", "HPV36L1":"AAA79442.1", "HPV36L2":"AAA79442.1", "HPV37E1":"AAA79445.1", "HPV37E2":"AAA79446.1", "HPV37E7":"AAA79444.1", "HPV37L1":"AAA79449.1", "HPV37L2":"AAA79449.1", "HPV38E1":"AAA79452.1", "HPV38E2":"AAA79453.1", "HPV38E7":"AAA79451.1", "HPV38L1":"AAA79456.1", "HPV38L2":"AAA79455.1", "HPV41E1":"NP_040287.1", "HPV41E2":"NP_040289.1", "HPV41E7":"NP_040286.1", "HPV41L1":"NP_040294.1", "HPV41L2":"NP_040293.1", "HPV47E1":"AAA46978.1", "HPV47E2":"AAA46979.1", "HPV47E7":"AAA46977.1", "HPV47L1":"AAA46982.1", "HPV47L2":"AAA46981.1", "HPV48E1":"NP_043418.1", "HPV48E2":"NP_043419.1", "HPV48E7":"NP_043417.1", "HPV48L1":"NP_043422.1", "HPV48L2":"NP_043421.1", "HPV49E1":"NP_041834.1", "HPV49E2":"NP_041835.1", "HPV49E7":"NP_041833.1", "HPV49L1":"NP_041837.1", "HPV49L2":"NP_041836.1", "HPV50E1":"NP_043425.1", "HPV50E2":"NP_043426.1", "HPV50E7":"NP_043424.1", "HPV50L1":"NP_043429.1", "HPV50L2":"NP_043428.1", "HPV60E1":"NP_043439.1", "HPV60E2":"NP_043440.1", "HPV60E7":"NP_043438.1", "HPV60L1":"NP_043443.1", "HPV60L2":"NP_043442.1", "HPV62E1":"AAR32248.1", "HPV62E2":"AAR32249.1", "HPV62E7":"AAR32247.1", "HPV62L1":"AAR32252.1", "HPV62L2":"AAR32251.1", "HPV63E1":"NP_040898.1", "HPV63E2":"NP_040899.1", "HPV63E7":"NP_040897.1", "HPV63L1":"NP_040902.1", "HPV63L2":"NP_040901.1", "HPV65E1":"CAA50173.1", "HPV65E2":"CAA50174.1", "HPV65E7":"CAA50172.1", "HPV65L1":"CAA50177.1", "HPV65L2":"CAA50176.1", "HPV67E1":"BAA66111.1", "HPV67E2":"BAA66112.1", "HPV67E7":"BAA66110.1", "HPV67L1":"BAA28859.1", "HPV67L2":"BAA66115.1", "HPV69E1":"BAA90729.1", "HPV69E2":"BAA90730.1", "HPV69E7":"BAA90728.1", "HPV69L1":"BAA90734.1", "HPV69L2":"BAA90733.1", "HPV74E1":"AAO15457.1", "HPV74E2":"AAO15458.1", "HPV74E7":"AAO15456.1", "HPV74L1":"AAO15462.1", "HPV74L2":"AAO15461.1", "HPV75E1":"CAA75451.1", "HPV75E2":"CAA75452.1", "HPV75E7":"CAA75450.1", "HPV75L1":"CAA75454.1", "HPV75L2":"CAA75453.1", "HPV76E1":"CAA75458.1", "HPV76E2":"CAA75459.1", "HPV76E7":"CAA75457.1", "HPV76L1":"CAA75461.1", "HPV76L2":"CAA75460.1", "HPV77E1":"CAA75465.1", "HPV77E2":"CAA75466.1", "HPV77E7":"CAA75464.1", "HPV77L1":"CAA75468.1", "HPV77L2":"CAA75467.1"}
proteinas2 = {}
for clave, valor in CargaBD.proteinas.items():
    proteinas2[valor] = clave


def ContruirArboles(elecion=0):
    global tree1
    datos_generados = []
    print("3: Se generaran arboles filogeneticos con los alineamientos grupales, espere por favor ...")
    archivo=archivos_MSA[elecion]
    #for archivo in archivos_MSA:
    # print("Creando arbol para:",archivo)
    archivo_salida = archivo.replace(".phylip", ".jpg")
    titulo = "Grupo: " + str(archivo[3:-13]).replace("_", " ") + " - Proteina: " + str(archivo[-13:-11])
    datos = titulo + " - salida: " + str(archivo_salida)
    datos_generados.append(datos)
    # construyo arbol Filogenetico
    # paso 1 leo arcchivos
    aln = AlignIO.read(archivo, 'phylip')
    # print(aln)
    # calculo de distancia con matriz default(dna_models, protein_models, models)
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(aln)
    # print(dm)
    # No se bien
    for count, value in enumerate(dm.names):
        if dicc_id_cepa.get(value):
            dm.names[count] = dicc_id_cepa[value]#[:-1]
        else:
            claves=dicc_id_cepa.keys()
            clave_sub= [x for x in claves if value in x]
            dm.names[count] = dicc_id_cepa[clave_sub[0]]#[:-1]
    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(aln)
    tree = constructor.upgma(dm)
    tree.ladderize()  # Imprime imagen
    # fig = plt.figure(figsize=(10, 20), dpi=100)
    # axes = fig.add_subplot(1, 1, 1)
    # Phylo.draw(tree, axes=axes, label_func=get_label)
    # fig=Phylo.draw(tree, label_func=lambda x: None, show_confidence=False, do_show=False)
    #fig = Phylo.draw(tree, do_show=True, branch_labels=proteinas2, show_confidence=True)
    fig = Phylo.draw(tree, do_show=False)
    plt.title(titulo)
    plt.ylabel("Cepas HPV")
    plt.xlabel("Distancias")
    plt.savefig(archivo_salida)
    plt.show()
    # print(tree.clade.clades)
    # for i in (tree.clade.clades):
    # print(i.clades)
    # for j in (i.clades):
    # print(j.name)
    # print()
    # print(tree.find_clades)
    # print(tree.clade.find_any)
    # print(tree.clade.find_clades)
    # print(tree.clade.find_elements)
    # print(dir(tree.clade))
    # Phylo.draw_ascii(tree) # Imprime por consola arbol
    print("Calculos terminados se han generado arboles filogenticos para:")
    for dato in datos_generados:
        print(dato)


def SubMenu1():
    print()
    print("### Analisis de resultados ###")
    print()
    print("Que desea realizar:")
    print("1- Graficar Proteina E1")
    print("2- Graficar Proteina E2")
    print("3- Graficar Proteina E7")
    print("4- Graficar Proteina L1")
    print("5- Graficar Proteina L2")
    print("6- Volver al menu anterior")
    print("7- Volver al menu principal")
    print("8- Salir")
    CargaBD.LeerInicio()
    opcion_principal = int(input("Ingrese una opción: "))
    if opcion_principal == 1:
        ContruirArboles(opcion_principal-1)
        SubMenu1()
    if opcion_principal == 2:
        ContruirArboles(opcion_principal-1)
        SubMenu1()
    if opcion_principal == 3:
        ContruirArboles(opcion_principal-1)
        SubMenu1()
    if opcion_principal == 4:
        ContruirArboles(opcion_principal-1)
        SubMenu1()
    if opcion_principal == 5:
        ContruirArboles(opcion_principal-1)
        SubMenu1()
    if opcion_principal == 6:
        Menu()
    if opcion_principal == 7:
        principal.Menu()
    if opcion_principal == 8:
        print("Gracias por utilizar BIOG5")
        sys.exit()

def Main():
    print("1: Obteniendo datos y agrupandolos por riesgo, espere por favor ...")
    ObtenerAgrupaciones()
    print("2: Se alineran genes-proteinas en grupos de riesgo, espere por favor ...")
    CorrerClustalOmega()
    LeerMSA()
    RelacionarIdCepa()
    SubMenu1()
    #ContruirArboles()
    Alineamientos.Limpiar(archivos_agrupados)  # ver si borrar xq requiere re alinear todo
    Alineamientos.Limpiar(archivos_agrupados_gen)
    Menu()
# Main()

def Menu():
    print()
    print("### Modulo: Alineamiento multiples ###")
    print()
    print("Que desea realizar:")
    print("1- Listar Proteinas almacenadas")
    print("2- Listar Genomas almacenadas")
    print("3- Realizar alineamientos multiples con Clustal Omega y analizar ")
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
    if opcion_principal == 4:
        principal.Menu()
    if opcion_principal == 5:
        print("Gracias por utilizar BIOG5")
        sys.exit()
# Menu()
