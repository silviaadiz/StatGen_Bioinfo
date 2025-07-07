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

1. **Paso 1**: Con PLINK, dividimos os datos por cromosomas (1-23) e convertémolos a formato VCF. É importante que os cromosomas se codifiquen co prefixo `chr` (ex. `chr1`, `chrX`) para que sexan compatibles coa nomenclatura do ficheiro FASTA e co servidor de imputación.

<aside>
👉

*NOTAS:*  En PLINK, a opción `--output-chr chrM` activa esta codificación (`chr1`, `chrX`, etc.). Ademais, PLINK 1.9 xera VCFs compatibles co estándar ≤4.2, que é o aceptado polos servidores de imputación.

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
    

1. **Paso 3:**  O servidor de imputación pode fallar cando hai demasiadas variantes con alelos desaxustados. Para evitar isto, usamos `bcftools norm` para eliminar as variantes cuxo alelo de referencia non coincide co xenoma. Os pasos concretos varían segundo a versión do xenoma:
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

Estes 23 arquivos se suben a imputar. O propio servidor de imputación fai un QC sobre os datos (https://topmedimpute.readthedocs.io/en/latest/pipeline/), e ten unha pipeline propia do chrX (revisar anterior link).

Unha vez completados os pasos anteriores, dispoñemos de 23 ficheiros `.vcf.gz` listos para a imputación. Deben cumprir os seguintes requisitos:

- **hg19**: Cromosomas sen prefixo `chr`, e `X` codificado como `23`.
- **hg38**: Cromosomas co prefixo `chr`, e `chrX` mantido como tal.

Estes ficheiros poden subirse ao servidor de imputación correspondente (por exemplo, TOPMed Imputation Server). O propio servidor executará un control de calidade previo, incluída unha pipeline específica para o cromosoma X. Pódese consultar máis información na documentación oficial:  [TOPMed Imputation Pipeline](https://topmedimpute.readthedocs.io/en/latest/pipeline/)

## Post-imputación

Unha vez completada a imputación, seguimos estes pasos para descargar, revisar e preparar os datos imputados para análises posteriores.

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

<aside>
👉

**Nota:** o alelo alternativo NON ten por que ser o menor. Podemos atopar algo de info aquí: https://genome.sph.umich.edu/wiki/Minimac3_Info_File.

**Nota:** Nós queremos traballar cos xenotipos codificados como 0/1/2, é dicir, os xenotipos ou campo GT. Se nalgún caso quixéramos recuperar as doses, poderíamos (https://www.cog-genomics.org/plink/2.0/input#vcf). 

</aside>

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

O nome dos SNPs orixinais está en formato rs, pero cambiámolo a chr:pos:ref:alt usando `--set-all-var-ids`.

<aside>
👉

ACLARACIÓNS IMPORTANTES

*Nomes SNPs e rs*

Os SNPs saen co rs, pero temos acceso ao chr:pos:ref:alt grazas ao arquivo info.gz, e o comando de plink `--set-all-var-ids '@:#:$r:$a'` permite construílo así. 

NOTA: Plink protesta se hai un indel moi grande, podemos arranxar isto cambiando o `--new-id-max-allele-len` (agora está en 100).

*Reincorporar FAM*

É crucial manter a coherencia nos datos dos individuos (FID, IID, sexo, fenotipo). Como PLINK asume que os individuos están na mesma orde, debemos asegurarnos de que os ficheiros `.fam` imputados coinciden cos orixinais xenotipados.

Se seguimos o código que hai no apartado Pre-imputación no debería haber problema, pois partimos dos mesmos individuos. Isto é especialmente importante en casos de tratar por separado algún dos cromosomas, como por exemplo o X: pode pasar que o QC sexa distinto, teñamos distintas versións do mesmo cromosoma, etc. e non partamos dos mesmos sample ID (errores pasados). 

</aside>

```bash
# Paso 2: Filtrado e renomeamento de SNPs dende VCF imputado
for i in {1..23}; do
    $STORE2/plink/plink2 \
        --vcf "${DIR_IMPUTACION}/chr${i}.dose.vcf.gz" \
        --set-all-var-ids '@:#:$r:$a' \
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

```

## Código completo

Aquí temos o código completo, pero no cesga están os arquivos “Script_pre_imputación.sh”, “Script_post_imputación.sh” e “[config.sh](http://config.sh)”, que son máis cómodos de usar porque só se modifica config.sh (salvo modificacións da pipeline).
