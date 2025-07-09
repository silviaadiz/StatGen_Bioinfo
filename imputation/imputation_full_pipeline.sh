#---------------------------------------
#---------- PRE IMPUTACIÓN -------------
#---------------------------------------

# Paso 0: módulos e variables

module load plink
module load gcccore/system samtools/1.9
DIR_XENOT="cambiar_a_directorio_xenotipados"
DIR_PRE_IMPUTACION="cambiar_a_directorio_output_pre_imputacion" 
# copiar o anterior se os queremos na mesma carpeta
NOME_XENOT="nome_arquivos_xenotipados"
GENOME_BUILD="hg19" 
# ou "hg38"
FASTA_hg19="/mnt/netapp1/Store_chumxrcg/REFERENCE_GENOMES/hg19.fa"
FASTA_hg38="/mnt/netapp1/Store_chumxrcg/REFERENCE_GENOMES/hg38.fa"
# Arquivo para renomear o código de cromosoma en hg19: chr1-->1; chrX-->23
RENAME_FILE="/mnt/netapp1/Store_chumxrcg/REFERENCE_GENOMES/RENAME_hg19_back.txt" 

# Paso 1: División en cromosomas e recodificación a VCF. Exportación co código de cromosoma correcto.

for i in {1..23}; do
    /mnt/netapp1/Store_chumxrcg/plink/plink \
        --bfile "${DIR_XENOT}/${NOME_XENOT}" \
        --chr "${i}" \
        --recode vcf-iid \
        --output-chr chrM \
        --out "${DIR_PRE_IMPUTACION}/${NOME_XENOT}_chr${i}"
done

# Paso 2: Compresión do VCF a VCF.gz e indexar

for i in {1..23}; do 
    bcftools view "${DIR_PRE_IMPUTACION}/${NOME_XENOT}_chr${i}.vcf" -Oz > "${DIR_PRE_IMPUTACION}/${NOME_XENOT}_chr${i}.vcf.gz"
    tabix -p vcf "${DIR_PRE_IMPUTACION}/${NOME_XENOT}_chr${i}.vcf.gz"
done

# Paso 3: Alineamento ao xenoma de referencia e cambio de código de chr en hg19

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


#---------------------------------------
#---------- POST IMPUTACIÓN ------------
#--------------------------------------- 

# Paso 0. Cargamos módulos e establecemos as variables.

module load plink
module load gcccore/system samtools/1.9
DIR_XENOT="cambiar_a_directorio_xenotipados"
DIR_PRE_IMPUTACION="cambiar_a_directorio_output_pre_imputacion" 
# copiar o anterior se os queremos na mesma carpeta
NOME_XENOT="nome_arquivos_xenotipados"
DIR_IMPUTACION="cambiar_path_para_descarga_imputados"
# Prefixo para os arquivos de saída
PREFIX_IMPUTADOS="cambiar_nome_para_imputados"
LINK_DESCARGA="copiar_link_descarga"
PASSWORD="copiar_password"

# Paso 1: Descarga dos arquivos dende o servidor
 
cd ${DIR_IMPUTACION}/
curl -sL ${LINK_DESCARGA} | bash
unzip -P ${PASSWORD} \*

# Paso 2: Filtrado e renomeamento de SNPs dende VCF imputado
for i in {1..23}; do
    $STORE2/plink/plink2 \
        --vcf "${DIR_IMPUTACION}/chr${i}.dose.vcf.gz" \
        --set-all-var-ids 'chr@:#:$r:$a' \
        --extract-if-info "R2>=0.8" \
        --new-id-max-allele-len 100 \
        --maf 0.01 \
        --allow-no-sex \
        --make-bed \
        --out "${DIR_IMPUTACION}/${PREFIX_IMPUTADOS}_chr${i}_post_qc_paso01"
done

# Paso 3: Reemprazo do arquivo FAM co orixinal pre-imputación
# NOTA: COMPROBAR QUE O CHR IMPUTADO CORRESPONDE AOS INDIVIDUOS XENOTIPADOS, é dicir, que sexan OS MESMOS INDIVIDUOS E NA MESMA ORDE

for i in {1..23}; do
    $STORE2/plink/plink \
        --bfile "${DIR_IMPUTACION}/${PREFIX_IMPUTADOS}_chr${i}_post_qc_paso01" \
        --fam "${DIR_XENOT}/${NOME_XENOT}.fam" \
        # NOTA: Se tiveramos os arquivos ${DIR_XENOT}/${NOME_XENOT}_chr${i}.fam/.bim/.bed (é dicir, crear no paso 1 pre-imputación os binarios por chr, ademais dos vcf)
        # poderíamos empalmarlle o fam do chr, é o mesmo. o comando sería substituír o de encima por:
        # --fam "${DIR_XENOT}/${NOME_XENOT}_chr${i}.fam"
        --allow-no-sex \
        --make-bed \
        --out "${DIR_IMPUTACION}/${PREFIX_IMPUTADOS}_chr${i}_imputado"
done
