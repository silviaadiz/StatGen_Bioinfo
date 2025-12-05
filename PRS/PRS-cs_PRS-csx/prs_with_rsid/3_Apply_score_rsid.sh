#!/bin/bash


# 0. Configuración de rutas e nomes

module load plink
PHENO="nome_feno"
TARGET_DATA_FOLDER="ruta_a_individuos_target"
INPUT_DATA_NAME="prefixo_nome_individuos_target"
SNPS_EXTRACT="arquivo_de_snps_a_extraer"

# SCORE_FILE é o nome que lle queremos poñer ao scorefile
SCORE_FILE="SNP.RSID.Scores...blablabla.txt"
OUTPUT_FOLDER="ruta_onde_gardar_output"
OUTPUT_SCORE_NAME="nome_output_prs"

mkdir "${OUTPUT_FOLDER}"


# 1. PASO DIFERENTE: Concatenar scores (non temos o paso de RSID to SNP)

for chr in {1..22}; do
  cat PRS_${PHENO}_pst_eff_a1_b0.5_phiauto_chr${chr}.txt >> "${SCORE_FILE}"
done

cut -f2 "${SCORE_FILE}" > ${SNPS_EXTRACT}


# 2. Extraemos os SNPs por cromosoma e eliminamos os rsid que teñan ".". Cambiamos a PLINK 2 por esta última función. quitamos  o "_imputado do final"
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

# 6. Eliminamos arquivos temporais
rm -f "${OUTPUT_FOLDER}"/temp.SNPS.PRScs_${INPUT_DATA_NAME}_chr*.* \
     "${OUTPUT_FOLDER}/temp.allchr.SNPS.PRScs_${INPUT_DATA_NAME}."* \
      "${OUTPUT_FOLDER}/list.files.merge.prscsx.${INPUT_DATA_NAME}.txt" \
