# 0. Configuración de rutas e nomes
PHENO="PTSD"
TARGET_DATA_FOLDER="/mnt/lustre/scratch/nlsas/home/usc/gb/sdd/EB25/PRS_GIULIA/00_imputados/"
INPUT_DATA_NAME="new_GIU_imputado08_mon"
SNPS_EXTRACT="extract_snps_prscs_${PHENO}_to_GIULIA.txt"
SCORE_FILE="SNP.RSID.Scores_prscs_${PHENO}_to_GIULIA.txt"
OUTPUT_FOLDER="/mnt/lustre/scratch/nlsas/home/usc/gb/sdd/EB25/PRS_GIULIA/01_PRS_psy/PRS_${PHENO}/"
OUTPUT_SCORE_NAME="${PHENO}_to_GIULIA"


module load plink

mkdir -p "${OUTPUT_FOLDER}"


# 1. Concatenar scores (non temos o paso de RSID to SNP)
OUT_NAME="OUT_NAME"

for chr in {1..22}; do
  cat ${OUT_NAME}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt >> ${SCORE_FILE}
done

cut -f2 ${SCORE_FILE} > ${SNPS_EXTRACT}

# 2. Extraemos os SNPs por cromosoma e eliminamos os rsid que teñan ".". Cambiamos a PLINK 2 por esta última función. 
for chr in {1..22}; do
    $LUSTRE/plink2 \
        --bfile "${TARGET_DATA_FOLDER}/${INPUT_DATA_NAME}_chr${chr}" \
        --rm-dup exclude-all \
        --extract "${SNPS_EXTRACT}" \
        --make-bed \
        --out "${OUTPUT_FOLDER}/temp.SNPS.PRScs_${INPUT_DATA_NAME}_chr${chr}" \
        --allow-no-sex
done

# 3. Creamos a lista de arquivos para mergear
find "${OUTPUT_FOLDER}" -name "temp.SNPS.PRScs_${INPUT_DATA_NAME}_chr*.bim" \
    | sort \
    | sed 's/\.bim$//' > "${OUTPUT_FOLDER}/list.files.merge.prscsx.${INPUT_DATA_NAME}.txt"

# 4. Facemos o merge dos cromosomas
plink \
    --bfile "${OUTPUT_FOLDER}/temp.SNPS.PRScs_${INPUT_DATA_NAME}_chr1" \
    --merge-list "${OUTPUT_FOLDER}/list.files.merge.prscsx.${INPUT_DATA_NAME}.txt" \
    --make-bed \
    --out "${OUTPUT_FOLDER}/temp.allchr.SNPS.PRScs_${INPUT_DATA_NAME}" \
    --allow-no-sex
    
# 5. Aplicamos o score final
plink \
    --bfile "${OUTPUT_FOLDER}/temp.allchr.SNPS.PRScs_${INPUT_DATA_NAME}" \
    --score "${SCORE_FILE}" 2 4 6 sum \
    --out "${OUTPUT_FOLDER}/PRS.PRScs.${OUTPUT_SCORE_NAME}" \
    --allow-no-sex
