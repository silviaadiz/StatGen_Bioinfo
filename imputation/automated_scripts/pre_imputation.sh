#!/bin/bash

#-----------------------------------------
# Paso 0.  Cargar configuración e módulos
#-----------------------------------------

module load plink
module load gcccore/system samtools/1.9
source config.sh

echo "Procesando con referencia $GENOME_BUILD"
echo "Arquivo de xenotipados orixinais está en $DIR_XENOT"
echo "Os arquivos xenotipados para subir a imputar estarán en $DIR_PRE_IMPUTACION e terán como prefixo $NOME_XENOT"


#----------------------------------------------------------------------------------------------------
# Paso 2: División en cromosomas e recodificación a VCF. Exportación co código de cromosoma correcto.
#----------------------------------------------------------------------------------------------------


for i in {1..23}; do
    /mnt/netapp1/Store_chumxrcg/plink/plink \
        --bfile "${DIR_XENOT}/${NOME_XENOT}" \
        --chr "${i}" \
        --recode vcf-iid \
        --output-chr chrM \
        --out "${DIR_PRE_IMPUTACION}/${NOME_XENOT}_chr${i}"
done

#----------------------------------------------
# Paso 3: Compresión do VCF a VCF.gz e indexar
#----------------------------------------------

for i in {1..23}; do 
    bcftools view "${DIR_PRE_IMPUTACION}/${NOME_XENOT}_chr${i}.vcf" -Oz > "${DIR_PRE_IMPUTACION}/${NOME_XENOT}_chr${i}.vcf.gz"
    tabix -p vcf "${DIR_PRE_IMPUTACION}/${NOME_XENOT}_chr${i}.vcf.gz"
done

#-------------------------------------------------------------------------------
# Paso 4: Alineamento ao xenoma de referencia e cambio de código de chr en hg19
#-------------------------------------------------------------------------------

for i in {1..23}; do
    if [ "$GENOME_BUILD" = "hg19" ]; then
    # Alineamento ao xenoma de referencia
        bcftools norm --check-ref s -f "$FASTA_hg19" "${DIR_PRE_IMPUTACION}/${NOME_XENOT}_chr${i}.vcf.gz" -o "${DIR_PRE_IMPUTACION}/chr${i}_flipped.vcf"
		# Reanotación dos SNPs (de chr1:pos:alelos a 1:pos:alelos)
        bcftools annotate --rename-chrs "$RENAME_FILE" "${DIR_PRE_IMPUTACION}/chr${i}_flipped.vcf" -Oz -o "${DIR_PRE_IMPUTACION}/chr${i}_flipped_preimputacion.vcf.gz"
    fi

    if [ "$GENOME_BUILD" = "hg38" ]; then
    # Alineamento ao xenoma de referencia
        bcftools norm --check-ref s -f "$FASTA_hg38" "${DIR_PRE_IMPUTACION}/${NOME_XENOT}_chr${i}.vcf.gz" -Oz -o "${DIR_PRE_IMPUTACION}/chr${i}_flipped_preimputacion.vcf.gz"
    fi
done
