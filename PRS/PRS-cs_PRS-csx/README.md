#!/bin/bash

module load plink

# 1. Configuración de rutas e nomes
TARGET_DATA_FOLDER="ruta/a/carpeta/con/binarios/imputados"
INPUT_DATA_NAME="nome_binarios_sin_chr"
SNPS_EXTRACT="archivo/snps/extraer/para/score"
SCORE_FILE="ruta/a/score/file"
OUTPUT_FOLDER="ruta/a/output"
OUTPUT_SCORE_NAME="nome_de_scores_aplicados"
# Asumindo que xa se 
OUT_NAME_PRSCS="OUT_NAME"

mkdir -p "${OUTPUT_FOLDER}"


# 2. Extraemos os SNPs por cromosoma
for chr in {1..22}; do
    plink \
        --bfile "${TARGET_DATA_FOLDER}/${INPUT_DATA_NAME}_chr${chr}_imputado" \
        --extract "${SNPS_EXTRACT}" \
        --make-bed \
        --out "${OUTPUT_FOLDER}/temp.SNPS.PRScs_${INPUT_DATA_NAME}_chr${chr}" \
        --allow-no-sex
done

# 3. Creamos a lista de arquivos para mergear
find "${OUTPUT_FOLDER}" -name "temp.SNPS.PRScs_${INPUT_DATA_NAME}_chr*.bim" \
    | sort \
    | sed 's/\.bim$//' > "${OUTPUT_FOLDER}/list.files.merge.prscs.${INPUT_DATA_NAME}.txt"

# 4. Facemos o merge dos cromosomas
plink \
    --bfile "${OUTPUT_FOLDER}/temp.SNPS.PRScs_${INPUT_DATA_NAME}_chr1" \
    --merge-list "${OUTPUT_FOLDER}/list.files.merge.prscs.${INPUT_DATA_NAME}.txt" \
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
      "${OUTPUT_FOLDER}/list.files.merge.prscs.${INPUT_DATA_NAME}.txt" \
