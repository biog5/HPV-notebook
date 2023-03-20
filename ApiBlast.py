#!/usr/bin/env python
# coding: utf-8

# # TP6 Grupo Nº 5
# Perez Ernesto Rafael rafaelperezctes@gmail.com Sofia Erdozain sofierdozain@gmail.com
# Created on Wed Nov  3 14:02:59 2021

# @author: Grupo Nº 5



import warnings
warnings.filterwarnings("ignore")
import requests
from io import StringIO
from Bio.Blast import NCBIXML

import requests
import time
import os

import requests
import time
import os


def run_blast(query_file, subject_file):
    base_url = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'

    # Leer el contenido de los archivos
    with open(query_file, 'r') as f:
        query_sequence = f.read()
    with open(subject_file, 'r') as f:
        subject_sequence = f.read()

    # Configurar los parámetros para la solicitud POST
    payload = {
        'CMD': 'Put',
        'PROGRAM': 'blastp',  # Asegúrate de usar 'blastn' para ADN
        'QUERY': query_sequence,
        'SUBJECTS': subject_sequence
    }

    # Enviar la solicitud POST
    response = requests.post(base_url, data=payload)
    response.raise_for_status()

    # Extraer el ID de la búsqueda (RID) de la respuesta
    result_page = response.text
    rid = None
    for line in result_page.splitlines():
        if line.startswith("RID = "):
            rid = line[6:]

    if rid is None:
        raise Exception("No se pudo obtener el RID del resultado")

    # Esperar a que el resultado esté listo
    payload = {
        'CMD': 'Get',
        'RID': rid
    }
    while True:
        time.sleep(10)  # Esperar 10 segundos antes de verificar nuevamente
        response = requests.post(base_url, data=payload)
        response.raise_for_status()
        result_page = response.text
        if "Status=WAITING" in result_page:
            continue
        elif "Status=READY" in result_page:
            break
        else:
            raise Exception("Error en el estado de la búsqueda")

    # Extraer la secuencia alineada de la respuesta y devolverla como una cadena de texto
    result_page = result_page.split('\n')
    alignment_started = False
    output_sequence = ''
    for line in result_page:
        if line.startswith('Sequences producing significant alignments:'):
            alignment_started = True
        elif alignment_started:
            if line.startswith('>'):
                output_sequence += '\n' + line + '\n'
            else:
                output_sequence += line.strip()

    return output_sequence


# En esta versión de la función, después de esperar a que el resultado esté listo, la función extrae la secuencia alineada de la respuesta del servidor de BLAST y la devuelve como una cadena de texto. La secuencia alineada se devuelve en formato FASTA, con una línea que comienza con el carácter '>' seguida del nombre de la secuencia y luego las líneas de la secuencia alineada.
# Para utilizar la función con tus archivos de entrada, puedes llamarla de la siguiente manera:

query_file = 'BD/HPV16E1.fasta'
subject_file = 'BD/low_riskE1.fasta'

#aligned_sequence = run_blast(query_file, subject_file)
#print(aligned_sequence)
import requests
import re
import time
from Bio.Blast import NCBIXML

# Definir el archivo de secuencia de consulta
query_file = 'BD/HPV16E1.fasta'


import requests
import re
import time
from Bio.Blast import NCBIXML

def blast_search(query_file):
    # Leer el archivo de secuencia de consulta
    with open(query_file, 'r') as f:
        query_sequence = f.read()

    # Construir el payload de la solicitud POST
    payload = {
        'CMD': 'Put',
        'PROGRAM': 'blastp',
        'DATABASE': 'nr',
        'QUERY': query_sequence,
        'EXPECT': 10.0
    }

    # Enviar la solicitud POST
    response = requests.post('https://blast.ncbi.nlm.nih.gov/Blast.cgi', data=payload)

    # Obtener el RID y el RTOE
    response_text = response.text
    rid_pattern = re.compile('RID = (.+)')
    rtoe_pattern = re.compile('RTOE = (.+)')
    rid_match = rid_pattern.search(response_text)
    rtoe_match = rtoe_pattern.search(response_text)

    if rid_match and rtoe_match:
        rid = rid_match.group(1)
        rtoe = rtoe_match.group(1)
        print(f'RID: {rid}')
        print(f'RTOE: {rtoe}')
    else:
        print('No se pudo obtener el RID o el RTOE de la página de resultados')

    # Esperar a que finalice la búsqueda
    while True:
        time.sleep(5)
        status_url = f'https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID={rid}'
        status_response = requests.get(status_url)
        status_text = status_response.text
        if 'Status=READY' in status_text:
            break

    # Obtener los resultados de la búsqueda
    result_url = f'https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_TYPE=XML&RID={rid}'
    result_response = requests.get(result_url)
    result_text = result_response.text
    with open('blast.xml', 'w') as f:
        f.write(result_text)

    # Procesar los resultados
    result_handle = open('blast.xml')
    blast_records = NCBIXML.parse(result_handle)
    blast_record = next(blast_records)
    best_hit = blast_record.alignments[0]
    print(f'Mejor coincidencia: {best_hit.hit_id}')

#blast_search('BD/HPV16E1.fasta')

import requests
import re
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.pairwise2 import align


def blast_searchB(query_file, subject_file):
    # Leer el archivo de secuencia de consulta
    with open(query_file, 'r') as f:
        query_sequence = f.read()

    # Leer el archivo de secuencias de sujeto
    with open(subject_file, 'r') as f:
        subject_sequence = f.read()

    # Construir el payload de la solicitud POST
    payload = {
        'CMD': 'Put',
        'PROGRAM': 'blastp',
        'QUERY': query_sequence,
        'SUBJECTS': subject_sequence,
        'EXPECT': 10.0
    }

    # Enviar la solicitud POST
    response = requests.post('https://blast.ncbi.nlm.nih.gov/Blast.cgi', data=payload)

    # Obtener el RID y el RTOE
    response_text = response.text
    rid_pattern = re.compile('RID = (.+)')
    rtoe_pattern = re.compile('RTOE = (.+)')
    rid_match = rid_pattern.search(response_text)
    rtoe_match = rtoe_pattern.search(response_text)

    if rid_match and rtoe_match:
        rid = rid_match.group(1)
        rtoe = rtoe_match.group(1)
        print(f'RID: {rid}')
        print(f'RTOE: {rtoe}')
    else:
        print('No se pudo obtener el RID o el RTOE de la página de resultados')

    # Esperar a que finalice la búsqueda
    while True:
        time.sleep(5)
        status_url = f'https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID={rid}'
        status_response = requests.get(status_url)
        status_text = status_response.text
        if 'Status=READY' in status_text:
            break

    # Obtener los resultados de la búsqueda
    result_url = f'https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_TYPE=XML&RID={rid}'
    result_response = requests.get(result_url)
    result_text = result_response.text
    with open('blast.xml', 'w') as f:
        f.write(result_text)

    # Procesar los resultados
    result_handle = open('blast.xml')
    blast_records = NCBIXML.parse(result_handle)
    blast_record = next(blast_records)
    best_hit = blast_record.alignments[0]
    print(f'Mejor coincidencia: {best_hit.hit_id}')

#blast_searchB('BD/HPV16E1.fasta', 'BD/low_riskE1.fasta')


def blast_searchC0():
    # Establecer las secuencias de consulta y sujeto
    query_sequence = 'GATCTCGAATCTCGGGTGCCAAGGAACTCCAGTCAC'
    subject_sequence = 'ATCTCGGGTGCCAAGGAACTCCAGTCACGATCTCGA'

    # Construir el payload de la solicitud POST
    payload = {
        'CMD': 'Put',
        'PROGRAM': 'blastn',
        'QUERY': query_sequence,
        'SUBJECTS': subject_sequence,
        'EXPECT': 10.0
    }

    # Enviar la solicitud POST
    response = requests.post('https://blast.ncbi.nlm.nih.gov/Blast.cgi', data=payload)

    # Verificar el código de estado de la respuesta
    if response.status_code != 200:
        print(f'Error: El código de estado de la respuesta fue {response.status_code}')
        return

    # Obtener el RID y el RTOE
    response_text = response.text
    rid_pattern = re.compile('RID = (.+)')
    rtoe_pattern = re.compile('RTOE = (.+)')
    rid_match = rid_pattern.search(response_text)
    rtoe_match = rtoe_pattern.search(response_text)

    # Verificar si se obtuvo el RID y el RTOE
    if not rid_match or not rtoe_match:
        # Verificar si hay algún mensaje de error en la respuesta
        error_pattern = re.compile('Error: (.+)')
        error_match = error_pattern.search(response_text)
        if error_match:
            print(f'Error: {error_match.group(1)}')
        else:
            # Mostrar la respuesta para ver si hay algún campo que pueda contener el RID
            print('Error: No se pudo obtener el RID o el RTOE de la página de resultados')
            print('Respuesta:')
            print(response_text)
        return

    rid = rid_match.group(1)
    rtoe = rtoe_match.group(1)
    print(f'RID: {rid}')
    print(f'RTOE: {rtoe}')

    # Esperar a que finalice la búsqueda
    while True:
        time.sleep(5)
        status_url = f'https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID={rid}'
        status_response = requests.get(status_url)
        status_text = status_response.text
        if 'Status=READY' in status_text:
            break

    # Obtener los resultados de la búsqueda
    result_url = f'https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_TYPE=XML&RID={rid}'
    result_response = requests.get(result_url)
    result_text = result_response.text
    with open('blast.xml', 'w') as f:
        f.write(result_text)

    # Procesar los resultados
    result_handle = open('blast.xml')
    blast_records = NCBIXML.parse(result_handle)
    blast_record = next(blast_records)
    best_hit = blast_record.alignments[0]
    print(f'Mejor coincidencia: {best_hit.hit_id}')

#blast_searchC0()
# Use BlastAlign.cgi to submit blast2seq searches pero no funca


from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def blast_searchD(query, subject):
    # Crear el registro de secuencia que contiene ambas secuencias
    combined_sequence = SeqRecord(seq=query + subject, id="combined_seq")
    sequences = [combined_sequence]

    # Realizar la búsqueda de Blast utilizando la librería BioPython con los parámetros deseados
    result_handle = NCBIWWW.qblast(
        program='blastn',
        database='nr',
        sequence=sequences,
        expect=10,
        word_size=11,
        alignments=1,
        hitlist_size=1,
        descriptions=0,
        entrez_query='(none)',
        megablast=True
    )

    # Procesar los resultados
    blast_records = NCBIXML.parse(result_handle)
    blast_record = next(blast_records)

    # Mostrar la información de la mejor coincidencia obtenida
    alignment = blast_record.alignments[0]
    hsps = alignment.hsps
    match = hsps[0]

    print(f"Alignment statistics for match #1:\nScore\tExpect\t\tIdentities\tGaps\tStrand\n{match.score}\t{match.expect}\t{match.identities}/{match.align_length} ({match.identities/match.align_length*100:.0f}%)\t\t{match.gaps}/{match.align_length} ({match.gaps/match.align_length*100:.0f}%)\t{match.frame}")
query = "GATCTCGAATCTCGGGTGCCAAGGAACTCCAGTCAC"
subject = "ATCTCGGGTGCCAAGGAACTCCAGTCACGATCTCGA"
#blast_searchD(query, subject)

