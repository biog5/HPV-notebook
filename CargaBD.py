#!/usr/bin/env python
# coding: utf-8

# # TP6 Grupo Nº 5
# Perez Ernesto Rafael rafaelperezctes@gmail.com Sofia Erdozain sofierdozain@gmail.com
# Created on Wed Nov  3 14:02:59 2021

# @author: Grupo Nº 5

from Bio import Entrez
from Bio import SeqIO
import Bio.SearchIO as bpio
import Bio.SeqIO as bsio
import matplotlib.pyplot as plt
import numpy as np
import Utils

Entrez.email = 'my_email@example.com'

"""
genomas={"HPV16":"NC_001526.4", "HPV18":"NC_001357.1", "HPV45":"X74479.1", "HPV31":"J04353.1", "HPV33":"M12732.1", "HPV35H":"X74477.1", "HPV39":"M62849.1", "HPV51":"M62877.1", "HPV52":"X74481.1", "HPV56":"X74483.1", "HPV58":"D90400.1", "HPV59":"X77858.1", "HPV68A":"DQ080079.1", "HPV68B":"FR751039.1", "HPV73":"X94165.1", "HPV82":"AB027021.1", "HPV23":"U31781.1", "HPV53":"NC_001593.1", "HPV66":"U31794.1", "HPV6B":"NC_001355.1", "HPV11":"M14119.1", "HPV40":"X74478.1", "HPV42":"M73236.1", "HPV43":"AJ620205.1", "HPV44":"U31788.1", "HPV54":"NC_001676.1", "HPV61":"NC_001694.1", "HPV70":"U21941.1", "HPV72":"KJ145795.1", "HPV81":"AJ620209.1", "HPV34":"NC_001587.1", "HPV57B":"U37537.1", "HPV57C":"AB361563.1", "HPV83":"AF151983.1", "HPV2":"MN605988.1", "HPV2A":"X55964.1","HPV3":"X74462.1","HPV5":"NC_001531.1", "HPV5":"NC_001457.1", "HPV7":"NC_001595.1", "HPV9":"NC_001596.1", "HPV10":"NC_001576.1", "HPV12":"X74466.1", "HPV13":"X62843.1", "HPV14B":"X74467.1", "HPV15":"X74468.1", "HPV17":"X74469.1", "HPV19":"X74470.1", "HPV20":"U31778.1", "HPV21":"U31779.1", "HPV22":"U31780.1", "HPV23":"U31781.1", "HPV24":"U31782.1", "HPV25":"X74471.1", "HPV27B":"AB211993.1", "HPV28":"U31783.1", "HPV29":"U31784.1", "HPV30":"NC_038889.1", "HPV32":"NC_001586.1", "HPV36":"U31785.1", "HPV37":"U31786.1", "HPV38":"U31787.1", "HPV41":"NC_001354.1", "HPV47":"M32305.1", "HPV48":"NC_001690.1", "HPV49":"NC_001591.1", "HPV50":"NC_001691.1", "HPV60":"NC_001693.1", "HPV62":"AY395706.1", "HPV63":"NC_001458.1", "HPV65":"X70829.1", "HPV67":"D21208.1", "HPV69":"AB027020.1", "HPV74":"AF436130.1", "HPV75":"Y15173.1", "HPV76":"Y15174.1", "HPV77":"Y15175.1"}
proteinas={"HPV16E1":"NP_041327.2", "HPV16E2":"NP_041328.1", "HPV16E7":"NP_041326.1", "HPV18E1":"NP_040312.1", "HPV18E2":"NP_040313.1", "HPV18E7":"NP_040311.1", "HPV45E1":"CAA52575.1", "HPV45E2":"CAA52576.1", "HPV45E7":"CAA52574.1", "HPV31E1":"AAA46952", "HPV31E2":"AAA46953.1", "HPV31E7":"AAA46951.1", "HPV33E1":"AAA46960.1", "HPV33E2":"AAA46961.1", "HPV33E7":"AAA46959.1", "HPV35HE1":"CAA52563.1", "HPV35HE2":"CAA52564.1", "HPV35HE7":"CAA52562.1", "HPV39E1":"AAA47052.1", "HPV39E2":"AAA47053.1", "HPV39E7":"AAA47051.1", "HPV52E1":"CAA52587.1", "HPV52E2":"CAA52588.1", "HPV52E7":"CAA52586.1", "HPV56E2":"CAA52598.1", "HPV56E7":"CAA52597.1", "HPV58E1":"BAA31847.1", "HPV58E2":"BAA31848.1", "HPV58E7":"BAA31846.1", "HPV59E1":"CAA54851", "HPV59E2":"CAA54852.1", "HPV59E7":"CAA54850.1", "HPV68AE1":"AAZ39493.1", "HPV68AE2":"AAZ39494.1", "HPV68AE7":"AAZ39492.1", "HPV68BE1":"CBY85069.1", "HPV68BE2":"CBY85070.1", "HPV68BE7":"CBY85068.1", "HPV73E1":"CAA63884.1", "HPV73E2":"CAA63885.1", "HPV73E7":"CAA63883.1", "HPV82E1":"BAA90737.1", "HPV82E2":"BAA90738.1", "HPV82E7":"BAA90736.1", "HPV23E1":"AAA79410.1", "HPV23E2":"AAA79411.1", "HPV23E7":"AAA79409.1", "HPV53E1A":"NP_597792.1","HPV53E1B":"NP_597793.1", "HPV53E2":"NP_041846.1", "HPV53E7":"NP_041845.1", "HPV66E1":"AAA79501.1", "HPV66E2":"AAA79502.1", "HPV66E7":"AAA79500.1", "HPV6BE1":"NP_040298.1", "HPV6BE2":"NP_040299.1", "HPV6BE7":"NP_040297.1","HPV11E1":"AAA46929.1", "HPV11E2":"AAA46930.1", "HPV11E7":"AAA46928.1", "HPV40E1":"CAA52569.1", "HPV40E2":"CAA52570.1", "HPV40E7":"CAA52568.1","HPV43E1":"CAF05785.1", "HPV43E2":"CAF05786.1", "HPV43E7":"CAF05784.1", "HPV44E1":"AAA79459.1", "HPV44E2":"AAA79460.1", "HPV44E7":"AAA79458.1", "HP54VE1":"NP_043290.1", "HPV54E2":"NP_043291.1", "HPV54E7":"NP_043289.1", "HPV61E1":"NP_043446.1", "HPV61E2":"NP_043447.1", "HPV61E7":"NP_043445.1", "HPV70E1":"AAC54852.1", "HPV70E2":"AAC54853.1", "HPV70E7":"AAC54851.1", "HPV72E1":"AHZ58190.1", "HPV72E2":"AHZ58191.1", "HPV72E7":"AHZ58194.1", "HPV81E1":"CAF05694.1", "HPV81E2":"CAF05695.1", "HPV81E7":"CAF05693.1", "HPV34E1":"NP_041809.1", "HPV34E2":"NP_041810.1", "HPV34E7":"NP_041808.1", "HPV57CE1":"BAF80482.1", "HPV57CE2":"BAF80483.1", "HPV57CE7":"BAF80481.1", "HPV83E1":"AAD38970.1", "HPV83E2":"AAD38971.1", "HPV83E7":"AAD38969.1", "HPV16L1":"NP_041332.2", "HPV16L2":"NP_041331.2", "HPV18L1":"NP_040317.1", "HPV18L2":"NP_040316.1", "HPV45L1":"CAA52578.1", "HPV45L2":"CAA52577.1", "HPV31L1":"AAA46956.1", "HPV31L2":"AAA46955.1", "HPV33L1":"AAA46964.1", "HPV33L2":"AAA46963.1", "HPV35HL1":"CAA52566.1", "HPV35HL2":"CAA52565.1", "HPV39L1":"AAA47056.1", "HPV39L2":"AAA47055.1", "HPV52L1":"CAA52590.1", "HPV52L2":"CAA52589.1", "HPV56L1":"CAA52600.1", "HPV56L2":"CAA52599.1", "HPV58L1":"BAA31851.1", "HPV58L2":"BAA31850.1", "HPV59L1":"CAA54856.1", "HPV59L2":"CAA54855.1", "HPV68AL1":"AAZ39498.1", "HPV68AL2":"AAZ39497.1", "HPV68BL1":"CBY85074.1", "HPV68BL2":"CBY85073.1", "HPV73L1":"CAA63887.1", "HPV73L2":"CAA63886.1", "HPV82L1":"BAA90742.1", "HPV82L2":"BAA90741.1", "HPV23L1":"AAA79414.1", "HPV23L2":"AAA79413.1", "HPV53L1":"NP_041848.1", "HPV53L2":"NP_041847.1", "HPV66L1":"AAA79505.1", "HPV66L2":"AAA79504.1", "HPV6BL1":"NP_040304.1", "HPV6BL2":"NP_040303.1", "HPV11L1":"AAA46935.1", "HPV11L2":"AAA46934.1", "HPV40L1":"CAA52572.1", "HPV40L2":"CAA52571.1","HPV42E1":"AAA47043.1", "HPV42E7":"AAA47042.1" ,"HPV42L1":"AAA47048.1", "HPV42L2":"AAA47047.1","HPV43L1":"CAF05790.1", "HPV43L2":"CAF05788.1", "HPV54L1":"NP_043294.1", "HPV54L2":"NP_043293.1", "HPV61L1":"NP_043450.1", "HPV61L2":"NP_043449.1", "HPV70L1":"AAC54857.1", "HPV70L2":"AAC54856.1", "HPV72L1":"AHZ58195.1", "HPV72L2":"AHZ58196.1", "HPV81L1":"CAF05698.1", "HPV81L2":"CAF05697.1","HPV34L1":"NP_041812.1", "HPV34L2":"NP_041811.1" ,"HPV57BL1":"AAC56600.1", "HPV57CL1":"BAF80486.1", "HPV57CL2":"BAF80485.1", "HPV83L1":"AAD38974.1", "HPV83L2":"AAD38973.1", "HPV2E1":"QLM04844.1", "HPV2E2":"QLM04841.1", "HPV2E7":"QLM04843.1", "HPV2L1":"QLM04846.1", "HPV2L2":"QLM04845.1", "HPV3E1":"CAA52471.1", "HPV3E2":"CAA52472.1", "HPV3L1":"CAA52474.1", "HPV3L2":"CAA52473.1","HPV5E1":"NP_041367.1", "HPV5E2":"NP_041368.1", "HPV5E7":"NP_041366.1", "HPV5L1":"NP_041372.1", "HPV5L2":"NP_041371.1", "HPV4E1":"NP_040891.1", "HPV4E2":"NP_040892.1", "HPV4E7":"NP_040890.1", "HPV4L1":"NP_040895.1", "HPV4L2":"NP_040894.1", "HPV7E1":"NP_041856.1", "HPV7E2":"NP_041857.1", "HPV7E7":"NP_041855.1", "HPV7L1":"NP_041859.1", "HPV7L2":"NP_041858.1", "HPV9E1":"NP_041863.1", "HPV9E2":"NP_041864.1", "HPV9E7":"NP_041862.1", "HPV9L1":"NP_041866.1", "HPV9L2":"NP_041865.1", "HPV10E1":"NP_041743.1", "HPV10E2":"NP_041744.1", "HPV10E7":"NP_041742.1", "HPV10L1":"NP_041746.1", "HPV10L2":"NP_041745.1", "HPV120E1":"CAA52498.1", "HPV12E2":"CAA52499.1", "HPV12E7":"CAA52497.1", "HPV12L1":"CAA52501.1", "HPV12L2":"CAA52500.1", "HPV13E1":"CAA44649.1", "HPV13E2":"CAA44650.1", "HPV13E7":"CAA44648.1", "HPV13L1":"CAA44654.1", "HPV13L2":"CAA44653.1", "HPV14BE1":"CAA52502.1", "HPV14BE2":"CAA52503.1", "HPV14BL1":"CAA52505.1", "HPV14BL2":"CAA52504.1", "HPV15E1":"CAA52508.1", "HPV15E2":"CAA52509.1", "HPV15E7":"CAA52507.1", "HPV15L1":"CAA52511.1", "HPV15L2":"CAA52510.1", "HPV17E1":"CAA52514.1", "HPV17E2":"CAA52515.1", "HPV17E7":"CAA52513.1", "HPV17L1":"CAA52517.1", "HPV17L2":"CAA52516.1", "HPV19E1":"CAA52520.1", "HPV19E2":"CAA52521.1", "HPV19E7":"CAA52519.1", "HPV19L1":"CAA52523.1", "HPV19L2":"CAA52522.1", "HPV20E1":"AAA79389.1", "HPV20E2":"AAA79390.1", "HPV20E7":"AAA79388.1", "HPV20L1":"AAA79393.1", "HPV20L2":"AAA79392.1", "HPV21E1":"AAA79396.1", "HPV21E2":"AAA79397.1", "HPV21E7":"AAA79395.1", "HPV21L1":"AAA79400.1", "HPV21L2":"AAA79399.1", "HPV22E1":"AAA79403.1", "HPV22E2":"AAA79404.1", "HPV22E7":"AAA79402.1", "HPV22L1":"AAA79407.1", "HPV22L2":"AAA79406.1", "HPV23E1":"AAA79410.1", "HPV23E2":"AAA79411.1", "HPV23E7":"AAA79409.1", "HPV23L1":"AAA79414.1", "HPV23L2":"AAA79413.1", "HPV24E1":"AAA79417.1", "HPV24E2":"AAA79418.1", "HPV24E7":"AAA79416.1", "HPV24L1":"AAA79421.1", "HPV24L2":"AAA79420.1", "HPV25E1":"CAA52526.1", "HPV25E2":"CAA52527.1", "HPV25E7":"CAA52525.1", "HPV25L1":"CAA52529.1", "HPV25L2":"CAA52528.1", "HPV27BE1":"BAE16265.1", "HPV27BE2":"BAE16266.1", "HPV27BE7":"BAE16264.1", "HPV27BL1":"BAE16269.1", "HPV27BL2":"BAE16268.1", "HPV28E1":"AAA79424.1", "HPV28E2":"AAA79425.1", "HPV28E7":"AAA79423.1", "HPV28L1":"AAA79428.1", "HPV28L2":"AAA79427.1", "HPV29E1":"AAA79431.1", "HPV29E2":"AAA79432.1", "HPV29E7":"AAA79430.1", "HPV29L1":"AAA79435.1", "HPV29L2":"AAA79434.1", "HPV30E1":"YP_009508156.1", "HPV30E2":"YP_009508157.1", "HPV30E7":"YP_009508155.1", "HPV30L1":"YP_009508159.1", "HPV30L2":"YP_009508158.1", "HPV32E1":"NP_041803.1", "HPV32E2":"NP_041804.1", "HPV32E7":"NP_041802.1", "HPV32L1":"NP_041806.1", "HPV32L2":"NP_041805.1", "HPV36E1":"AAA79438.1", "HPV36E2":"AAA79439.1", "HPV36E7":"AAA79437.1", "HPV36L1":"AAA79442.1", "HPV36L2":"AAA79442.1", "HPV37E1":"AAA79445.1", "HPV37E2":"AAA79446.1", "HPV37E7":"AAA79444.1", "HPV37L1":"AAA79449.1", "HPV37L2":"AAA79449.1", "HPV38E1":"AAA79452.1", "HPV38E2":"AAA79453.1", "HPV38E7":"AAA79451.1", "HPV38L1":"AAA79456.1", "HPV38L2":"AAA79455.1", "HPV41E1":"NP_040287.1", "HPV41E2":"NP_040289.1", "HPV41E7":"NP_040286.1", "HPV41L1":"NP_040294.1", "HPV41L2":"NP_040293.1", "HPV47E1":"AAA46978.1", "HPV47E2":"AAA46979.1", "HPV47E7":"AAA46977.1", "HPV47L1":"AAA46982.1", "HPV47L2":"AAA46981.1", "HPV48E1":"NP_043418.1", "HPV48E2":"NP_043419.1", "HPV48E7":"NP_043417.1", "HPV48L1":"NP_043422.1", "HPV48L2":"NP_043421.1", "HPV49E1":"NP_041834.1", "HPV49E2":"NP_041835.1", "HPV49E7":"NP_041833.1", "HPV49L1":"NP_041837.1", "HPV49L2":"NP_041836.1", "HPV50E1":"NP_043425.1", "HPV50E2":"NP_043426.1", "HPV50E7":"NP_043424.1", "HPV50L1":"NP_043429.1", "HPV50L2":"NP_043428.1", "HPV60E1":"NP_043439.1", "HPV60E2":"NP_043440.1", "HPV60E7":"NP_043438.1", "HPV60L1":"NP_043443.1", "HPV60L2":"NP_043442.1", "HPV62E1":"AAR32248.1", "HPV62E2":"AAR32249.1", "HPV62E7":"AAR32247.1", "HPV62L1":"AAR32252.1", "HPV62L2":"AAR32251.1", "HPV63E1":"NP_040898.1", "HPV63E2":"NP_040899.1", "HPV63E7":"NP_040897.1", "HPV63L1":"NP_040902.1", "HPV63L2":"NP_040901.1", "HPV65E1":"CAA50173.1", "HPV65E2":"CAA50174.1", "HPV65E7":"CAA50172.1", "HPV65L1":"CAA50177.1", "HPV65L2":"CAA50176.1", "HPV67E1":"BAA66111.1", "HPV67E2":"BAA66112.1", "HPV67E7":"BAA66110.1", "HPV67L1":"BAA28859.1", "HPV67L2":"BAA66115.1", "HPV69E1":"BAA90729.1", "HPV69E2":"BAA90730.1", "HPV69E7":"BAA90728.1", "HPV69L1":"BAA90734.1", "HPV69L2":"BAA90733.1", "HPV74E1":"AAO15457.1", "HPV74E2":"AAO15458.1", "HPV74E7":"AAO15456.1", "HPV74L1":"AAO15462.1", "HPV74L2":"AAO15461.1", "HPV75E1":"CAA75451.1", "HPV75E2":"CAA75452.1", "HPV75E7":"CAA75450.1", "HPV75L1":"CAA75454.1", "HPV75L2":"CAA75453.1", "HPV76E1":"CAA75458.1", "HPV76E2":"CAA75459.1", "HPV76E7":"CAA75457.1", "HPV76L1":"CAA75461.1", "HPV76L2":"CAA75460.1", "HPV77E1":"CAA75465.1", "HPV77E2":"CAA75466.1", "HPV77E7":"CAA75464.1", "HPV77L1":"CAA75468.1", "HPV77L2":"CAA75467.1"}

plantar_warts = ["HPV1", "HPV2", "HPV4", "HPV63"]
common_warts = ["HPV2", "HPV1", "HPV7", "HPV4", "HPV26", "HPV27", "HPV29", "HPV41", "HPV57", "HPV65", "HPV77", "HPV1", "HPV3", "HPV4", "HPV10", "HPV28"]
flat_warts = ["HPV3", "HPV10", "HPV26", "HPV27", "HPV28", "HPV38", "HPV41", "HPV49", "HPV75", "HPV76"]
other_cutaneous_lesions = ["HPV6", "HPV11", "HPV16", "HPV30", "HPV33", "HPV36", "HPV37", "HPV38", "HPV41", "HPV48", "HPV60", "HPV72", "HPV73"]
epidermodysplasia_verruciformis = ["HPV2", "HPV3", "HPV10", "HPV5", "HPV8", "HPV9", "HPV12", "HPV14", "HPV15", "HPV17", "HPV19", "HPV20", "HPV21", "HPV22", "HPV23", "HPV24", "HPV25", "HPV36", "HPV37", "HPV38", "HPV47", "HPV50"]
recurrent_respiratory_papillomatosis = ["HPV6", "HPV11"]
focal_epithelial_hyperplasia_of_heck = ["HPV13", "HPV32"]
conjunctival_papillomas_carcinomas = ["HPV6", "HPV11", "HPV16"]
condyloma_acuminata = ["HPV6", "HPV11", "HPV30", "HPV42", "HPV43", "HPV45", "HPV51", "HPV54", "HPV55", "HPV70"]
#Cervical intraepithelial neoplasia

unspecified = ["HPV30", "HPV34", "HPV39", "HPV40", "HPV53", "HPV57", "HPV59", "HPV61", "HPV62", "HPV64", "HPV66", "HPV67", "HPV68", "HPV69"]
low_risk = ["HPV6", "HPV11", "HPV16", "HPV18", "HPV31", "HPV33", "HPV35", "HPV42", "HPV43", "HPV44", "HPV45", "HPV51", "HPV52", "HPV74"]
high_risk = ["HPV16", "HPV18", "HPV6", "HPV11", "HPV31", "HPV34", "HPV33", "HPV35", "HPV39", "HPV42", "HPV44", "HPV45", "HPV51", "HPV52", "HPV56", "HPV58", "HPV66"]
cervical_carcinoma = ["HPV16", "HPV18", "HPV31", "HPV45", "HPV33", "HPV35", "HPV39", "HPV51", "HPV52", "HPV56", "HPV58", "HPV66", "HPV68", "HPV70"]
"""
import csv

genomas = {}
unspecified = None
high_risk = None
low_risk = None
proteinas = {}
proteinas_dictReader = {}
genomas_dictReader = {}
archivos_agrupados_gen = ['BD/E1.fasta', 'BD/E2.fasta', 'BD/E7.fasta', 'BD/L1.fasta', 'BD/L2.fasta']
archivos_agrupados = ['BD/high_riskE1.fasta', 'BD/high_riskE2.fasta', 'BD/high_riskE7.fasta', 'BD/high_riskL1.fasta',
                      'BD/high_riskL2.fasta', 'BD/high_riskE1.fasta', 'BD/high_riskE2.fasta', 'BD/high_riskE7.fasta',
                      'BD/high_riskL1.fasta', 'BD/high_riskL2.fasta', 'BD/high_riskE1.fasta', 'BD/high_riskE2.fasta',
                      'BD/high_riskE7.fasta', 'BD/high_riskL1.fasta', 'BD/high_riskL2.fasta', 'BD/high_riskE1.fasta',
                      'BD/high_riskE2.fasta', 'BD/high_riskE7.fasta', 'BD/high_riskL1.fasta', 'BD/high_riskL2.fasta',
                      'BD/high_riskE1.fasta', 'BD/high_riskE2.fasta', 'BD/high_riskE7.fasta', 'BD/high_riskL1.fasta',
                      'BD/high_riskL2.fasta', 'BD/high_riskE1.fasta', 'BD/high_riskE2.fasta', 'BD/high_riskE7.fasta',
                      'BD/high_riskL1.fasta', 'BD/high_riskL2.fasta', 'BD/high_riskE1.fasta', 'BD/high_riskE2.fasta',
                      'BD/high_riskE7.fasta', 'BD/high_riskL1.fasta', 'BD/high_riskL2.fasta', 'BD/high_riskE1.fasta',
                      'BD/high_riskE2.fasta', 'BD/high_riskE7.fasta', 'BD/high_riskL1.fasta', 'BD/high_riskL2.fasta',
                      'BD/high_riskE2.fasta', 'BD/high_riskE7.fasta', 'BD/high_riskL1.fasta', 'BD/high_riskL2.fasta',
                      'BD/high_riskE1.fasta', 'BD/high_riskE2.fasta', 'BD/high_riskE7.fasta', 'BD/high_riskL1.fasta',
                      'BD/high_riskL2.fasta', 'BD/low_riskE1.fasta', 'BD/low_riskE2.fasta', 'BD/low_riskE7.fasta',
                      'BD/low_riskL1.fasta', 'BD/low_riskL2.fasta', 'BD/low_riskE1.fasta', 'BD/low_riskE2.fasta',
                      'BD/low_riskE7.fasta', 'BD/low_riskL1.fasta', 'BD/low_riskL2.fasta', 'BD/low_riskE1.fasta',
                      'BD/low_riskE7.fasta', 'BD/low_riskL1.fasta', 'BD/low_riskL2.fasta', 'BD/low_riskE1.fasta',
                      'BD/low_riskE2.fasta', 'BD/low_riskE7.fasta', 'BD/low_riskL1.fasta', 'BD/low_riskL2.fasta',
                      'BD/low_riskE1.fasta', 'BD/low_riskE2.fasta', 'BD/low_riskE7.fasta', 'BD/low_riskL1.fasta',
                      'BD/low_riskL2.fasta', 'BD/low_riskE2.fasta', 'BD/low_riskE7.fasta', 'BD/low_riskL1.fasta',
                      'BD/low_riskL2.fasta', 'BD/low_riskE1.fasta', 'BD/low_riskE2.fasta', 'BD/low_riskE7.fasta',
                      'BD/low_riskL1.fasta', 'BD/low_riskL2.fasta', 'BD/unspecified_riskE1.fasta',
                      'BD/unspecified_riskE2.fasta', 'BD/unspecified_riskE7.fasta', 'BD/unspecified_riskL1.fasta',
                      'BD/unspecified_riskL2.fasta', 'BD/unspecified_riskE1.fasta', 'BD/unspecified_riskE2.fasta',
                      'BD/unspecified_riskE7.fasta', 'BD/unspecified_riskL1.fasta', 'BD/unspecified_riskL2.fasta',
                      'BD/unspecified_riskE1.fasta', 'BD/unspecified_riskE2.fasta', 'BD/unspecified_riskE7.fasta',
                      'BD/unspecified_riskL1.fasta', 'BD/unspecified_riskL2.fasta', 'BD/unspecified_riskE1.fasta',
                      'BD/unspecified_riskE2.fasta', 'BD/unspecified_riskE7.fasta', 'BD/unspecified_riskL1.fasta',
                      'BD/unspecified_riskL2.fasta', 'BD/unspecified_riskE1.fasta', 'BD/unspecified_riskE2.fasta',
                      'BD/unspecified_riskE7.fasta', 'BD/unspecified_riskL1.fasta', 'BD/unspecified_riskL2.fasta']
conservados = {1: ["Alto Riesgo", "E1", "high_riskE1_MSA"], 2: ["Alto Riesgo", "E2", "high_riskE2_MSA"],
               3: ["Alto Riesgo", "E7", "high_riskE7_MSA"], 4: ["Alto Riesgo", "L1", "high_riskL1_MSA"],
               5: ["Alto Riesgo", "L2", "high_riskL2_MSA"], 6: ["Bajo Riesgo", "E1", "low_riskE1_MSA"],
               7: ["Bajo Riesgo", "E2", "low_riskE2_MSA"], 8: ["Bajo Riesgo", "E7", "low_riskE7_MSA"],
               9: ["Bajo Riesgo", "L1", "low_riskL1_MSA"], 10: ["Bajo Riesgo", "L2", "low_riskL2_MSA"],
               11: ["Noespecificado Riesgo", "E1", "unspecified_riskE1_MSA"],
               12: ["Noespecificado Riesgo", "E2", "unspecified_riskE2_MSA"],
               13: ["Noespecificado Riesgo", "E7", "unspecified_riskE7_MSA"],
               14: ["Noespecificado Riesgo", "L1", "unspecified_riskL1_MSA"],
               15: ["Noespecificado Riesgo", "L2", "unspecified_riskL2_MSA"]}


def LeerListaProteinas():
    bd = 'proteinas.csv'
    archivo = "BD/" + bd
    global proteinas, proteinas_dictReader
    handle = open(archivo)
    reader = list(csv.DictReader(handle))
    proteinas_dictReader = reader
    # print("Leyendo lista de proteina, espere por favor...")
    for row in reader:
        key = row['#Genoma']
        valor = row['IdProtein']
        proteinas[key] = valor
    handle.close()


def LeerListaGenomas():
    # print("Ingrese en nombre completo del archivo (Debe estar en la carpeta BD/)")
    # import CargaBD
    bd = 'genomas.csv'
    archivo = "BD/" + bd
    handle = open(archivo)
    global genomas, genomas_dictReader
    reader = list(csv.DictReader(handle))
    genomas_dictReader = reader
    # print("Leyendo lista de genomas espere por favor...")
    for row in reader:
        # print(row['#Genoma'],row['IdProtein'])
        key = row['#Genoma']
        valor = row['IdGenoma']
        genomas[key] = valor
    handle.close()


def LeerListaClasificacion():
    bd = 'clasificacion.csv'
    archivo = "BD/" + bd
    global unspecified, high_risk, low_risk
    handle = open(archivo)
    reader = list(csv.reader(handle))
    # print("Leyendo lista de clasificacion espere por favor...")
    lista = {}
    for row in reader:
        lista[row[0]] = row[1:]
    handle.close()
    unspecified = lista["unspecified"]
    high_risk = lista["high_risk"]
    low_risk = lista["low_risk"]


def LeerInicio():
    # Mejorar
    if unspecified == None or high_risk == None or low_risk == None:
        LeerListaProteinas()
        LeerListaGenomas()
        LeerListaClasificacion()
        # print("aqui5")


# LeerInicio()

"""
#Aclaraciones
#Correr 2 veces si falla conexion
# Revordar Tener creada la carpeta "BD"

"""


############################# 1: Descargas

# Permite descargar genomas
def DescargarGenomas(output_file, records_to_download):
    for record_id in records_to_download:  # Deje el for para descargar todo junto
        handle = Entrez.efetch(db='nucleotide', id=record_id, rettype='gb')
        seqRecord = SeqIO.read(handle, format='gb')
        handle.close()
        output_file.write(seqRecord.format('fasta'))


# Permite descargar genomas
def DescargarProteinas(output_file, records_to_download):
    for record_id in records_to_download:
        handle = Entrez.efetch(db="protein", id=record_id, rettype="gb", retmode="text")
        seqRecord = SeqIO.read(handle, format='gb')
        handle.close()
        output_file.write(seqRecord.format('fasta'))


# Lee lista de genomas y lo envia a descargar en archivos fasta por separados
def GenomasSeparados():
    for v, k in genomas.items():
        archivo = 'BD/' + str(v) + '.fasta'
        if Utils.Existe(archivo) == False:
            output_file = open(archivo, "w")
            DescargarGenomas(output_file, [str(k)])
            print("Descargando: ", v, k)


# Lee lista de proteinas y lo envia a descargar en archivos fasta por separados
def ProteinasSeparadas():
    for v, k in proteinas.items():
        archivo = 'BD/' + str(v) + '.fasta'
        if Utils.Existe(archivo) == False:
            output_file = open(archivo, "w")
            DescargarProteinas(output_file, [str(k)])
            print("Descargando: ", v, k)


# Lee lista de genomas y lo envia a descargar todas juntos en un archivo fasta
def GenomasJuntos():
    lista_aux = []
    for v, k in genomas.items():
        lista_aux.append(k)
    archivo = 'BD/all_genomas.fasta'
    if Utils.Existe(archivo) == False:
        output_file = open(archivo, "w")
        print("Descargando: ", lista_aux)
        print("........")
        DescargarGenomas(output_file, lista_aux)


# Lee lista de proteinas y lo envia a descargar todas junto en un archivo fasta
def ProteinasJuntas():
    lista_aux = []
    archivo = 'BD/all_proteinas.fasta'
    if Utils.Existe(archivo) == False:
        for v, k in proteinas.items():
            lista_aux.append(k)
        output_file = open(archivo, "w")
        print("Descargando: ", lista_aux)
        print("........")
        DescargarProteinas(output_file, lista_aux)


def Main():
    GenomasSeparados()
    ProteinasSeparadas()
    GenomasJuntos()
    ProteinasJuntas()
# Main()