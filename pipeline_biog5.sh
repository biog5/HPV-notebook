#!/bin/bash
echo “----------------PIPELINE biog5--------------”

echo Obeniendo semillas para cada E1-alto riesgo
mafft high_riskE1.fasta > high_riskE1_seed_msa.fasta 
echo Obeniendo modelo para cada E1-alto riesgo
hmmbuild high_riskE1_modelo.hmm high_riskE1_seed_msa.fasta
echo Obeniendo secuencias de E1-bajo riesgo que concidan con cada modelo E1-alto riesgo
hmmsearch high_riskE1_modelo.hmm low_riskE1.fasta > high_riskE1_vs_low_riskE1_.out
echo Obeniendo secuencias de E1-no determinado riesgo que concidan con cada modelo E1-alto riesgo
hmmsearch high_riskE1_modelo.hmm unspecified_riskE1.fasta > high_riskE1_vs_unspecified_riskE1_.out


echo Obeniendo semillas para cada E2-alto riesgo
mafft high_riskE2.fasta > high_riskE2_seed_msa.fasta
echo Obeniendo modelo para cada E2-alto riesgo
hmmbuild high_riskE2_modelo.hmm high_riskE2_seed_msa.fasta
echo Obeniendo secuencias de E2-bajo riesgo que concidan con cada modelo E2-alto riesgo
hmmsearch high_riskE2_modelo.hmm low_riskE2.fasta > high_riskE2_vs_low_riskE2_.out
echo Obeniendo secuencias de E2-no determinado riesgo que concidan con cada modelo E2-alto riesgo
hmmsearch high_riskE2_modelo.hmm unspecified_riskE2.fasta > high_riskE2_vs_unspecified_riskE2_.out

echo Obeniendo semillas para cada E7-alto riesgo
mafft high_riskE7.fasta > high_riskE7_seed_msa.fasta
echo Obeniendo modelo para cada E7-alto riesgo
hmmbuild high_riskE7_modelo.hmm high_riskE7_seed_msa.fasta
echo Obeniendo secuencias de E7-bajo riesgo que concidan con cada modelo E7-alto riesgo
hmmsearch high_riskE7_modelo.hmm low_riskE7.fasta > high_riskE7_vs_low_riskE7_.out
echo Obeniendo secuencias de E7-no determinado riesgo que concidan con cada modelo E7-alto riesgo
hmmsearch high_riskE7_modelo.hmm unspecified_riskE7.fasta > high_riskE7_vs_unspecified_riskE7_.out

echo Obeniendo semillas para cada L1-alto riesgo
mafft high_riskL1.fasta > high_riskL1_seed_msa.fasta
echo Obeniendo modelo para cada L1-alto riesgo
hmmbuild high_riskL1_modelo.hmm high_riskL1_seed_msa.fasta
echo Obeniendo secuencias de L1-bajo riesgo que concidan con cada modelo L1-alto riesgo
hmmsearch high_riskL1_modelo.hmm low_riskL1.fasta > high_riskL1_vs_low_riskL1_.out
echo Obeniendo secuencias de L1-no determinado riesgo que concidan con cada modelo L1-alto riesgo
hmmsearch high_riskL1_modelo.hmm unspecified_riskL1.fasta > high_riskL1_vs_unspecified_riskL1_.out

echo Obeniendo semillas para cada L2-alto riesgo
mafft high_riskL2.fasta > high_riskL2_seed_msa.fasta
echo Obeniendo modelo para cada L2-alto riesgo
hmmbuild high_riskL2_modelo.hmm high_riskL2_seed_msa.fasta
echo Obeniendo secuencias de L2-bajo riesgo que concidan con cada modelo L2-alto riesgo
hmmsearch high_riskL2_modelo.hmm low_riskL2.fasta > high_riskL2_vs_low_riskL2_.out
echo Obeniendo secuencias de L2-no determinado riesgo que concidan con cada modelo L2-alto riesgo
hmmsearch high_riskL2_modelo.hmm unspecified_riskL2.fasta > high_riskL2_vs_unspecified_riskL2_.out



