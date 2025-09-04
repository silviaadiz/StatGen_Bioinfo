# Pipeline completa imputación

## Pre-imputación

Partimos dos datos post-QC en formato binario de PLINK (cromosomas 1 a 23), que conteñen unicamente os SNPs xenotipados. Para poder subilos ao servidor de imputación, é necesario convertelos a formato VCF.

1. **Paso 0:** cargamos as librerías en bash (samtools, PLINK 1.9) e establecemos os directorios e a versión do xenoma.

```bash
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
```

1. **Paso 1**: Con PLINK, dividimos os datos por cromosomas (1-23) e convertémolos a formato VCF. 

Cando traballamos coa versión hg38, o número de cromosoma debe ir codificado como `chr` (ex. `chr1`, `chrX`) para o servidor de imputación, mentres que coa hg19 non.
Para alinear co ficheiro FASTA no Paso 3, por outro lado, a codificación dos cromosomas ten que levar o prefixo para que sexan compatibles coa nomenclatura do ficheiro. Entón o que facemos na pipeline é incorporar o "chr" directamente e: (A) eliminalo nun paso posterior se a versión do xenoma é hg19, ou (B). mantelo se é hg38.

<aside>


**NOTAS:**  En PLINK, a opción `--output-chr chrM` activa esta codificación (`chr1`, `chrX`, etc.). Ademais, PLINK 1.9 xera VCFs compatibles co estándar ≤4.2, que é o aceptado polos servidores de imputación.

</aside>

```bash
# Paso 1: División en cromosomas e recodificación a VCF. Exportación co código de cromosoma correcto.

for i in {1..23}; do
    /mnt/netapp1/Store_chumxrcg/plink/plink \
        --bfile "${DIR_XENOT}/${NOME_XENOT}" \
        --chr "${i}" \
        --recode vcf-iid \
        --output-chr chrM \
        --out "${DIR_PRE_IMPUTACION}/${NOME_XENOT}_chr${i}"
done
```

1. **Paso 2:** Comprimimos os ficheiros VCF a formato `.vcf.gz` e xeramos os seus índices con `tabix`:

    
    ```bash
    # Paso 2: Compresión do VCF a VCF.gz e indexar
    
    for i in {1..23}; do 
        bcftools view "${DIR_PRE_IMPUTACION}/${NOME_XENOT}_chr${i}.vcf" -Oz > "${DIR_PRE_IMPUTACION}/${NOME_XENOT}_chr${i}.vcf.gz"
        tabix -p vcf "${DIR_PRE_IMPUTACION}/${NOME_XENOT}_chr${i}.vcf.gz"
    done
    ```
    

1. **Paso 3:**  O servidor de imputación pode fallar cando hai demasiadas variantes con alelos desalineados. Para evitar isto, usamos `bcftools norm --check-ref s`u, que comproba se os alelo de referencia do VCF e do fasta están alineados e corrixe os que non. 

Ademais, hai que corrixir a codificación dos cromosomas para o servidor:

- **hg19**: Tras a validación, é necesario eliminar o prefixo `chr` e converter `chrX` a `23`.
- **hg38**: Mantemos o formato con prefixo `chr`.

O parámetro `--check-ref s` elimina as variantes incompatibles co alelo de referencia.

```bash
# Paso 3: Alineamento ao xenoma de referencia e cambio de código de chr en hg19

# Arquivo para renomear o código de cromosoma en hg19: chr1-->1; chrX-->23
RENAME_FILE="/mnt/netapp1/Store_chumxrcg/REFERENCE_GENOMES/RENAME_hg19_back.txt" 

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

```

## Imputación

Unha vez completados os pasos anteriores, dispoñemos de 23 ficheiros `.vcf.gz` listos para a imputación. Deben cumprir os seguintes requisitos:

- **hg19**: Cromosomas sen prefixo `chr`, e `X` codificado como `23`.
- **hg38**: Cromosomas co prefixo `chr`, e `chrX` mantido como tal.

Estes ficheiros poden subirse ao servidor de imputación correspondente (por exemplo, TOPMed Imputation Server). O propio servidor executará un control de calidade previo, incluída unha pipeline específica para o cromosoma X. Pódese consultar máis información na documentación oficial:  [TOPMed Imputation Pipeline](https://topmedimpute.readthedocs.io/en/latest/pipeline/)

## Post-imputación

Unha vez completada a imputación, seguimos os seguintes pasos para descargar, revisar e preparar os datos imputados para análises posteriores.

1. **Paso 0:** Cargamos os módulos e establecemos as variables.
2. **Paso 1:** Descargamos**.** Tras a finalización da imputación, recibimos un correo co *link* e contrasinal de descarga. Cada cromosoma está dividido en dous ficheiros:
- `chrN.dose.vcf.gz`: Xenotipos imputados en formato VCF, con campos como **GT** (genotipo), **DS** (dosificación), e **GP** (probabilidades xenotípicas).
- `chrN.info.gz`: Contén métricas por variante, como R2, MAF, AF, etc.

```bash
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
LINK_DESCARGA="copiar_link"
PASSWORD="copiar_password"

# Paso 1: Descarga dos arquivos dende o servidor
 
cd ${DIR_IMPUTACION}
curl -sL ${LINK_DESCARGA} | bash
unzip -P ${PASSWORD} \*
```

Exemplo de encabezado en `dose.vcf.gz`:

```markdown
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##filedate=20250415
##phasing=full
##contig=<ID=chr12>
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated Alternate Allele Frequency">
##INFO=<ID=MAF,Number=1,Type=Float,Description="Estimated Minor Allele Frequency">
##INFO=<ID=AVG_CS,Number=1,Type=Float,Description="Average Call Score">
##INFO=<ID=R2,Number=1,Type=Float,Description="Estimated Imputation Accuracy (R-square)">
##INFO=<ID=ER2,Number=1,Type=Float,Description="Empirical (Leave-One-Out) R-square (available only for genotyped variants)">
##INFO=<ID=IMPUTED,Number=0,Type=Flag,Description="Marker was imputed">
##INFO=<ID=TYPED,Number=0,Type=Flag,Description="Marker was genotyped">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=A,Type=Float,Description="Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]">
##FORMAT=<ID=GP,Number=G,Type=Float,Description="Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1">
##FORMAT=<ID=HDS,Number=.,Type=Float,Description="Estimated Haploid Alternate Allele Dosage">
```


**Nota:** o alelo alternativo NON ten por que ser o menor. Podemos atopar algo de info aquí: https://genome.sph.umich.edu/wiki/Minimac3_Info_File.

**Nota:** Nós queremos traballar cos xenotipos codificados como 0/1/2, é dicir, os xenotipos ou campo GT. Se nalgún caso quixéramos recuperar as doses, poderíamos (https://www.cog-genomics.org/plink/2.0/input#vcf). 



Exemplo de encabezado en `info.gz`:

```markdown
##fileformat=VCFv4.2
##filedate=20250415
##source=Minimac v4.1.6
##phasing=full
##contig=<ID=chr22>
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated Alternate Allele Frequency">
##INFO=<ID=MAF,Number=1,Type=Float,Description="Estimated Minor Allele Frequency">
##INFO=<ID=AVG_CS,Number=1,Type=Float,Description="Average Call Score">
##INFO=<ID=R2,Number=1,Type=Float,Description="Estimated Imputation Accuracy (R-square)">
##INFO=<ID=ER2,Number=1,Type=Float,Description="Empirical (Leave-One-Out) R-square (available only for genotyped variants)">
##INFO=<ID=IMPUTED,Number=0,Type=Flag,Description="Marker was imputed">
##INFO=<ID=TYPED,Number=0,Type=Flag,Description="Marker was genotyped">
```

1. **Paso 2:** Convertimos os ficheiros VCF a binario PLINK (`.bed/.bim/.fam`), renomeamos os SNPs co formato `chr:pos:ref:alt`, e aplicamos os seguintes filtros:
- Imputación de calidade: `R2 >= 0.8`
- Filtro de frecuencia: `MAF >= 0.01`
    
    **Paso 2.b:** Incorporamos o fam correspondente ao cromosoma xenotipado para recuperar toda a info de FID, IID, sexo e phenotype. 
    


**ACLARACIÓNS IMPORTANTES**

*Nomes SNPs e rs*

Os SNPs saen co rs, pero algúns deles non teñen un identificador asignado e levan “.” como nome. Temos acceso ao chr:pos:ref:alt grazas ao arquivo info.gz, que xunto ao comando de plink `--set-all-var-ids '@:#:$r:$a'` nos permiten construír un formato `chr:pos:ref:alt`. 

NOTA: PLINK protesta se hai un indel moi grande, podemos arranxar isto modificando  o `--new-id-max-allele-len` (agora está en 100).

*Reincorporar FAM*

É crucial manter a coherencia nos datos dos individuos (FID, IID, sexo, fenotipo). Como PLINK asume que os individuos están na mesma orde, debemos asegurarnos de que os ficheiros `.fam` imputados coinciden cos orixinais xenotipados.

Se seguimos o código que hai no apartado *Pre-imputación* no debería haber problema, pois partimos dos mesmos individuos que imputamos. Isto é especialmente importante en caso de tratar por separado algún dos cromosomas, como por exemplo o X: pode pasar que o QC sexa distinto, teñamos distintas versións do mesmo cromosoma, etc. e non partamos dos mesmos sample ID (errores pasados). 

*Dosages en DS* 

GT son os xenotipos discretos ou *hard-calls* (xenotipo máis probable). **Sen embargo, nalgún caso a imputación ten unha probabilidade ao redor do 50% de cada alelo, polo que sería interesante exportalo como “DS”, que son as probabilidades dos xenotipos. Ao convertilo en binario de PLINK, convírtense en GT, polo que habería que exportalo en pgen. Por defecto, PLINK pon a missing os xenotipos que non teñen GT, pero podemos forzalo a saltarse esas  variantes con `--vcf-require-gt`

---

</aside>

```bash
# Paso 2: Filtrado e renomeamento de SNPs dende VCF imputado

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
        --make-pgen dosage=DS \
        --fam "${DIR_XENOT}/${NOME_XENOT}.fam" \
        --out "${DIR_IMPUTACION}/${PREFIX_IMPUTADOS}_chr${i}_imputado"
done

# NOTA: COMPROBAR QUE O CHR IMPUTADO CORRESPONDE AOS INDIVIDUOS XENOTIPADOS, é dicir, que sexan OS MESMOS INDIVIDUOS E NA MESMA ORDE
```

Se quixéramos exportar dosages:

```r
for i in {1..23}; do
    $STORE2/plink/plink2 \
        --vcf dosage=DS "${DIR_IMPUTACION}/chr${i}.dose.vcf.gz" \
        --set-all-var-ids 'chr@:#:$r:$a' \
        --extract-if-info "R2>=0.8" \
        --new-id-max-allele-len 100 \
        --maf 0.01 \
        --allow-no-sex \
        --make-pgen \
        --fam "${DIR_XENOT}/${NOME_XENOT}.fam" \
        --out "${DIR_IMPUTACION}/${PREFIX_IMPUTADOS}_chr${i}_imputado"
done

```

## Código completo

Aquí temos o código completo, pero no cesga están os arquivos “Script_pre_imputación.sh”, “Script_post_imputación.sh” e ”config.sh”, que son máis cómodos de usar porque só se modifica config.sh (salvo modificacións da pipeline).

```bash
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
 
cd ${DIR_IMPUTACION}
curl -sL ${LINK_DESCARGA} | bash
unzip -P ${PASSWORD} \*

# Paso 2: Filtrado e renomeamento de SNPs dende VCF imputado

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
        --make-pgen dosage=DS \
        --fam "${DIR_XENOT}/${NOME_XENOT}.fam" \
        --out "${DIR_IMPUTACION}/${PREFIX_IMPUTADOS}_chr${i}_imputado"
done

# NOTA: COMPROBAR QUE O CHR IMPUTADO CORRESPONDE AOS INDIVIDUOS XENOTIPADOS, é dicir, que sexan OS MESMOS INDIVIDUOS E NA MESMA ORDE

```


## Código completo

[Aquí temos o código completo](https://github.com/silviaadiz/scriptsGMX/blob/main/imputation/imputation_full_pipeline.sh), pero podemos atopar un [modelo automatizado](https://github.com/silviaadiz/scriptsGMX/tree/main/imputation/automated_scripts), cos arquivos “Script_pre_imputación.sh”, “Script_post_imputación.sh” e “config.sh”, que son máis cómodos de usar porque só se toca o script config.sh (salvo modificacións da pipeline).
