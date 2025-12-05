#!/bin/bash


cat /ruta/a/scores/_pst_eff_a1_b0.5_phiauto_chr* > scores_PRScs_SA_Doch_EUR.prscs_to_EBeur.and.GIULIA_allchr.txt

Rscript --vanilla convert_rsid_a_SNP.R

cut -f2 -d$'\t' "SNP.RSID.Scores_prscs_SA.Doch.EUR_to_EBeur.and.GIULIA.txt" > extract_snps_prscs_SA.Doch.EUR_to_EBeur.and.GIULIA.txt
