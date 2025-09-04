#!/bin/bash

#-------------------------------
# Paso 0.  Cargar configuración 
#------------------------------- 

source config.sh


echo "Os arquivos imputados estarán en $DIR_IMPUTACION e terán como prefixo $PREFIX_IMPUTADOS"

#-----------------------------------------------
# Paso 1: Descarga dos arquivos dende o servidor
#-----------------------------------------------
 
cd ${DIR_IMPUTACION}/
curl -sL ${LINK_DESCARGA} | bash

unzip -P ${PASSWORD} \*

#------------------------------------------------------------
# Paso 2a: Filtrado e renomeamento de SNPs dende VCF imputado
# Paso 2b: Reemprazo do arquivo FAM co orixinal pre-imputación
# NOTA: COMPROBAR QUE O CHR IMPUTADO CORRESPONDE AOS INDIVIDUOS XENOTIPADOS, é dicir, que sexan OS MESMOS INDIVIDUOS E NA MESMA ORDE
#------------------------------------------------------------

for i in {1..23}; do
    if [ "$i" -eq 23 ]; then
        VCF_FILE="chrX.dose.vcf.gz"
    else
        VCF_FILE="chr${i}.dose.vcf.gz"
    fi

    $STORE2/plink/plink2 \
        --vcf "${DIR_IMPUTACION}/${VCF_FILE}" \
        --set-all-var-ids 'chr@:#:$r:$a' \
        --extract-if-info "R2>=0.8" \
        --new-id-max-allele-len 100 \
        --maf 0.01 \
        --allow-no-sex \
        --make-bed \
        --fam "${DIR_XENOT}/${NOME_XENOT}.fam" \
        --out "${DIR_IMPUTACION}/${PREFIX_IMPUTADOS}_chr${i}_imputado"
done
