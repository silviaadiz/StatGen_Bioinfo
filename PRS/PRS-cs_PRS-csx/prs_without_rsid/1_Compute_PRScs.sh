# Ruta aos xenomas de referencia


REF_DIR="/mnt/lustre/scratch/nlsas/home/usc/gb/sdd/PRSs_IMPORTANTE/prscs_prscsx/ref_pops/ldblk_1kg_eur"

# Arquivo BIM creado anteriormente sen ".bim"
BIM_PREFIX="/mnt/lustre/scratch/nlsas/home/usc/gb/sdd/EB25/PRSs_novembro25_cohortes_combi/bim.rs.GIULIA.EBeur.allchr"

# Summary statistics formateados
SUMSTATS="/mnt/lustre/scratch/nlsas/home/usc/gb/sdd/EB25/PRS_Docherty/00_archivos_input_PRS/Sumstats/Prepared_PRScs_EUR_SuicideAttempt_Docherty_ISGC.MVP_310725.txt"

# N do GWAS (tamaño de mostra)
N_GWAS="815178"

# CHR que estamos usando (DESTA FORMA PARA PARALELIZAR)
CHR=${SLURM_ARRAY_TASK_ID}

# Directorio de output
OUT_DIR="/mnt/lustre/scratch/nlsas/home/usc/gb/sdd/EB25/PRSs_novembro25_cohortes_combi/PRS_Docherty/"

module load miniconda3/4.11.0


python /mnt/lustre/scratch/nlsas/home/usc/gb/sdd/PRSs_IMPORTANTE/prscs_prscsx/PRScs/PRScs.py \
  --ref_dir=${REF_DIR} \
  --bim_prefix=${BIM_PREFIX} \
  --sst_file=${SUMSTATS} \
  --n_gwas=${N_GWAS} \
  --chrom=${CHR} \
  --out_dir=${OUT_DIR}
