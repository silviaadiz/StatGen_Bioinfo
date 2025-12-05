#
# 1. Configuración de rutas e nomes
source config_applyscore.sh
module load plink

mkdir -p "${OUTPUT_FOLDER}"

# 2. Extraemos os SNPs por cromosoma
for chr in {1..22}; do
    plink \
        --bfile "${TARGET_DATA_FOLDER}/${INPUT_DATA_NAME}_chr${chr}" \
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
    --score "${SCORE_FILE}" 2 5 7 sum \
    --out "${OUTPUT_FOLDER}/PRS.PRScs.${OUTPUT_SCORE_NAME}" \
    --allow-no-sex

# 6. Eliminamos arquivos temporais
rm -f "${OUTPUT_FOLDER}"/temp.SNPS.PRScs_${INPUT_DATA_NAME}_chr*.* \
      "${OUTPUT_FOLDER}/temp.allchr.SNPS.PRScs_${INPUT_DATA_NAME}."* \
      "${OUTPUT_FOLDER}/list.files.merge.prscsx.${INPUT_DATA_NAME}.txt" \
