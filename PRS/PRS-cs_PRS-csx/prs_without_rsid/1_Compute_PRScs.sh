#!/bin/bash

# Ruta aos xenomas de referencia
REF_DIR="/ruta_prscs/ref_pops/ldblk_1kg_eur"

# Arquivo BIM creado anteriormente sen ".bim"
BIM_PREFIX="bim.rs.GIULIA.EBeur.allchr"

# Summary statistics formateados
SUMSTATS="/ruta/sumstats.txt"

# N do GWAS (tamaño de mostra)
N_GWAS="815178"

# CHR que estamos usando (DESTA FORMA PARA PARALELIZAR)
CHR=${SLURM_ARRAY_TASK_ID}

# Directorio de output
OUT_DIR="/ruta/"
PHENO_NAME="blabla"
module load miniconda3/4.11.0


python /mnt/lustre/scratch/nlsas/home/usc/gb/sdd/PRSs_IMPORTANTE/prscs_prscsx/PRScs/PRScs.py \
  --ref_dir=${REF_DIR} \
  --bim_prefix=${BIM_PREFIX} \
  --sst_file=${SUMSTATS} \
  --n_gwas=${N_GWAS} \
  --chrom=${CHR} \
  --out_dir="${OUT_DIR}/PRS_${PHENO_NAME}"
