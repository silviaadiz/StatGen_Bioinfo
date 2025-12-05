#!/bin/bash

PHENO_NAME="blabla"
SCOREFILE_OUTPUT="nome_para_scorefile"
TARGET_NAME="nome_cohorte_individuos_target_para_output"


cat /ruta/a/scores/PRS_${PHENO_NAME}_pst_eff_a1_b0.5_phiauto_chr* > "${SCOREFILE_OUTPUT}"

Rscript --vanilla convertir_rsid_a_SNP.R \
--equivalence_file "/ruta/a/equiv.file.rs.snp/.txt" --score "${SCOREFILE_OUTPUT}" --pheno_name ${PHENO_NAME} --target_name ${TARGET_NAME}
  
cut -f2 -d$'\t' "SNP.RSID.Scores_prscs_${PHENO_NAME}_to_${TARGET_NAME}.txt" > "extract_snps_prscs_${PHENO_NAME}_to_${TARGET_NAME}.txt"
