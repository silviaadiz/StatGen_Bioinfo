# **PRS-CS e PRS-CSx**

**PRS-CS** e **PRS-CSx** son uns métodos baiesianos para computar **PRS**.

No seguinte link está a wiki:

[🔗 https://github.com/getian107/PRScs](https://github.com/getian107/PRScs?tab=readme-ov-file#readme)

[🔗](https://github.com/getian107/PRScs?tab=readme-ov-file#readme)https://github.com/getian107/PRScsx

## **Descarga de datos:**

Para empregar este software precisamos:

1. **Clonar** o repositorio para descargar os scripts.
2. **Descargar** os paneis de referencia.
    
    En caso de querer usar **PRS-CSx** (multi-poboación), debemos descargar **TODOS** os paneis de referencia correspondentes a ditas poboacións.
    
    Pódense descargar os do **UKBB** ou os de **1000G**.
    
    O panel de referencia debe corresponderse ao panel de referencia do GWAS de discovery.
    
    > NOTA: Os paneis de referencia NON levan o chr23 e só se inclúen os SNPs de HapMap.
    > 
3. **Summary statistics** de interese.

### PRS-CS vs PRS-CSx

- **PRS-CS**
    - Usa un único arquivo de *summary statistics*.
    - Pode ser dunha soa poboación ou multi-población.
    - Se é multi-población, débese elixir o panel **LD** da poboación máis numerosa (ver https://www.sciencedirect.com/science/article/pii/S2666979X23002227) ou crear un novo panel.
- **PRS-CSx**
    - Toma *sumstats* de poboacións individuais e acopla os efectos.
    - Permite:
        - “Meta-analizar” as poboacións.
        - Integrar os *scores* como combinación lineal. Tras xerar os scores, aplicaríanse a unha cohorte *testing* para estimar os pesos para cada PRS.

## Paso 1. Preparar sumstats y bim.

### 01. Preparar sumstats

Unha vez temos descargados os *summary statistics*, debemos revisar os campos (beta, SE, pvalue, nome do SNP, etc). Este software SÓ quere as columnas: `SNP, A1, A2, SE , BETA/OR` ou `SNP, A1, A2, P, BETA/OR`.

`A1` é o alelo de efecto.

PRS-CS(x) traballa con *rsid*, polo que tanto o .bim dos datos diana como os sumstats deben ter unha columna *rsid*. 

En caso de que existan tanto a columna SNP como *rsid* nos sumstats do GWAS de descubrimento (o habitual), **empregaríamos o script `0_get_rsid.R`:

```r
sumstats<-read.table("ruta/a/sumstats/sumstats.txt",header=T)
sumstats2<-sumstats[c("rsid","A1","A2","BETA","SE")]
# En caso de que teñan outros nomes de columnas, seleccionaríanse
	# as que foran e cambiaríanse os nomes a SNP, A1, A2, BETA, P.

names(sumstats2)[1]<-"SNP" # DEBE CHAMARSE SNP
write.table(sumstats2,"ruta/a/carpeta/prs/sumstats_preparado.txt",quote=F,row.names=F)
```

<aside>
💡

Se non se teñen os *rsid* nos *sumstats* nin temos forma de obtelos doutro arquivo, podemos anotalos con ANNOVAR. Pode que os sumstats necesiten máis modificacións, isto vai depender de como estean.  Se é para o PRS-CSx, debemos preparalos todos. 

</aside>

### 02. Preparar bim.

Necesitamos un `.bim` correspondente aos datos diana para indicarlle ao software que SNPs temos dispoñibles nos nosos datos. En principio, se os datos están imputados en topmed, temos acceso aos rsid correspondentes ao [SNP.ID](http://SNP.ID) nos INFO. 

En principio, se os datos están imputados en **TOPMed**, temos acceso aos *rsid* e aos [SNP.ID](http://SNP.ID) (`chr:pos:ref:alt`) nos .INFO.

***Nota***: Realmente este `.bim` non se usará en PLINK, polo que dá igual eliminar SNPs sen tocar o `.bed`, pero para PRS-cs debemos manter o formato clásico. 

1. **En caso de que os datos veñan de TOPMed e se usaran os meus scripts de post-imputación:** copiamos a función `anotar_bim_con_rsid()` a un script que se chame `anotar_bim_con_rsid_topmed.R` e cargámola como está a continuación no script onde traballemos (así podemos reutilizala para varios proxectos). Outra opción sería copiar a función directamente no script de traballo, pero é máis cómodo tela nunha carpeta externa e chamala cando a queiramos usar.

***NOTA**: S*e non se usaron os meus scripts post-imputación, hai que ter coidado con esta función, especialmente ao ler o arquivo “.bim”:`ficheiro_bim <- paste0(bim_path_prefix,"/", i, "_imputado.bim"` Esta liña de código espera que os imputados leven o sufixo “*imputado.bim*”, e dará problemas se non se chaman igual. PRS-cs e PRS-csx non teñen en conta o chrX.

```r
anotar_bim_con_rsid_topmed <- function(chr = 1:22,
                                 bim_path_prefix = "/ruta/a/binarios/imputados/prefixo",
                                 info_path = "/ruta/a/info/",
																 out_path = "/ruta/a/carpeta/prs/",
                                 out_prefix = "nome_output") {
	 if (!dir.exists(out_path)) {
    dir.create(out_path, recursive = TRUE)
	}
  
  for (i in chr) {
    
    # 1. Ler o ficheiro BIM do cromosoma i
    ficheiro_bim <- paste0(bim_path_prefix, "chr", i, "_imputado.bim")
    bim <- read.table(ficheiro_bim)
    names(bim) <- c("CHR", "snp", "cm", "POS", "A1", "A2")
    
    # 2. Ler o ficheiro INFO correspondente
    
    ficheiro_info <- paste0(info_path, "/chr", i, ".info.gz")
    
    info <- read.table(ficheiro_info)
    names(info) <- c("CHR", "POS", "RSID", "REF", "ALT", "QUAL", "FILTER", "INFO")

    # 3. Crear identificador SNP no formato chr:pos:ref:alt
    info$SNP <- paste0(info$CHR, ":", info$POS, ":", info$REF, ":", info$ALT)
    
    # 4. Facer merge co BIM a través do identificador SNP
    bim_rs <- merge(bim, info[c("RSID", "SNP")], by.x = "snp", by.y = "SNP", sort = FALSE)
    
    # 5. Filtrar rexistros con rsid válido
    snps_validos <- subset(bim_rs, !is.na(RSID) & RSID != ".")
    
    
    # 6. Gardar o novo ficheiro BIM con rsIDs
    exportar_bim <- snps_validos[c("CHR", "RSID", "cm", "POS", "A1", "A2")]
		output_completo <- paste0(out_path, "/bim_rs_", out_prefix, "_chr", i, ".bim")
    write.table(exportar_bim, output_completo, quote = F, row.names = F, col.names = F)
		
   # 7. Gardamos as conversións rsid-SNP porque as necesitaremos despois
		exportar_conv <- snps_validos[c("RSID", "snp")]
		output_completo_2 <- paste0(out_path, "/equivalencia.rs.snp.", out_prefix, "_chr", i, ".txt")
		write.table(exportar_conv,  output_completo_2, quote = F, row.names = F, col.names = F)
}
}

```

Unha vez teñamos o script coa función nalgún lado, podemos cargalo así:

```r
source("ruta/a/script/anotar_bim_con_rsid_topmed.R")

# Executamos a función cos parámetros
anotar_bim_con_rsid_topmed(bim_path_prefix = "ruta/a/binarios/imputados/con/prefixo", 
										info_path = "ruta/a/info",
										out_path = "/ruta/a/prs", 
										out_prefix = "poñer_nome_output")
# NOTA: a ruta ao bim debe levar tamén o prefixo do bim (é dicir, todo excepto "_chr"). 
#       Exemplo:
#				"/ruta/EB25_EUR_chr22_imputado"
```

1. **En caso de que teñamos os *rsid* nos *sumstats*, poderíamos facer algo así:**

```r
# Abrimos sumstats
rsidst<-read.table("sumstats.txt")

# Procesar cromosomas do 1 ao 22
for (i in 1:22) {
  
  # Ler arquivo BIM do cromosoma i
  bim <- read.table(paste0("/ruta/a/binarios/imputados/prefixo", i, ".bim"))
  # Esta liña debe modificarse acorde ao nome que usamos para os cromosomas
    
  # Engadimos nomes ao bim.
  names(bim) <- c("CHR", "snp", "cm", "POS", "A1", "A2")

  # Merge da columna "snp" do BIM coa columna SNP dos sumstats
  bim_rs <- merge(bim, rsidst[c("rs","SNP")], by.x = "snp", by.y = "SNP",sort=F)
  
  # Filtramos as columnas para quedarnos coas necesarias no BIM
  exportar <- subset(bim_rs, !is.na(rs) & rs != ".", select = c("CHR", "rs", "cm", "POS", "A1", "A2"))

  # Gardamos novo BIM anotado con rsIDs
  write.table(exportar, paste0("bim_rs_chr", i, ".bim"), quote = FALSE, row.names = FALSE, col.names = FALSE)}
```

Por último, concatenamos os `.bim` en bash.

```bash
OUT_PREFIX="out_prefix_de_anotar_bim_con_rsid_topmed"

for chr in {1..22}; do
  cat bim_rs_${OUT_PREFIX}_chr${chr}.bim >> bim_rs_${OUT_PREFIX}_allchr.bim
  cat equivalencia.rs.snp.${OUT_PREFIX}_chr${chr}.txt >> equivalencia.rs.snp.${OUT_PREFIX}_allchr.txt
done
```

## Paso 2. Executamos PRS-CS

### 01. PRS-CS cunha única poboación ou meta-GWAS cross-ancestry

Con PRS-Cs construímos PRS dunha única poboación ou multi-poboación cunha meta-análise (pero cun só panel de LD de referencia). Para cada cromosoma, creamos un script .sh:

```bash

**#!/bin/bash**
#SBATCH -p short -t 3:00:00 
#SBATCH --mem=16G 
#SBATCH -n1
#SBATCH --mail-type=end
#SBATCH --mail-user=silviadiz03@gmail.com
#SBATCH --output=./logs/%A_%a_%x_%j.out
#SBATCH --cpus-per-task=1
#SBATCH --array=1-22

# Ruta aos xenomas de referencia
REF_DIR="/mnt/lustre/scratch/nlsas/home/usc/gb/sdd/PRSs_IMPORTANTE/prscsx/ref_pops/ldblk_1kg_eur"

# Arquivo BIM creado anteriormente sen ".bim"
BIM_PREFIX="/ruta/bim_rs_allchr"

# Summary statistics formateados
SUMSTATS="/ruta/sumstats/formateados"

# N do GWAS (tamaño de mostra)
N_GWAS="1000"

# CHR que estamos usando
CHR="1"

# Directorio de output
OUT_DIR="/ruta/output/"
OUT_NAME="out_name"

module load miniconda3/4.11.0
source activate prscsx

python /mnt/lustre/scratch/nlsas/home/usc/gb/sdd/PRSs_IMPORTANTE/prscs_prscsx/PRScs/PRScs.py \
  --ref_dir=${REF_DIR} \
  --bim_prefix=${BIM_PREFIX} \
  --sst_file=${SUMSTATS} \
  --n_gwas=${N_GWAS} \
  --chrom=${CHR} \
  --out_dir="${OUT_DIR}/${OUT_NAME}"

```

Ademais destas opcións, ten outros parámetros como os `a` e `b` dos *prior gamma,* ou o parámetro ***phi*** (*shrinkage*). Se non se especifica, o propio software estima *phi*, pero pódese facer un grid search con varios valores e probar cal resulta en maior poder preditivo. Wang e colaboradores (2023) determinaron que o rendemento da búsqueda automática de *phi* non dista significativamente do *grid search*. 

Tamén se poden pedir os *posterior effect sizes* para os xenotipos estandarizados, indicar as iteracións de *burn-in* para o algoritmo MCMC…etc.

📄 [Ver README oficial de PRS-CS](https://github.com/getian107/PRScs/blob/master/README.md)

### 02. PRS-CSx

PRS-CSx estima os efectos dos SNPs acoplándoos entre poboacións. Precisa sumstats individualizados por poboación. 
📄 [Tutorial PRS-CSx](https://github.com/getian107/PRScsx)

Ten unha opción `-meta`que meta-analiza os *posterior effect sizes* das diferentes poboacións. O ideal, como se comentou arriba, sería obter un PRS por cada poboación (con esta opción desactivada), aplicala a unha cohorte para estimar os beta asociados aos PRS, e aplicalos a unhas segunda (training) e terceira (test) cohortes.

En xeral, pode que haxa cromosomas que tarden moito (máis de 6h para o 1 ou 2). Para revisar fácilmente se todos acaban, engadimos ao final unhas liñas de código.

```bash
#!/bin/bash
#SBATCH -p short -t 6:00:00 
#SBATCH --mem=16G 
#SBATCH -n1
#SBATCH --mail-type=end
#SBATCH --mail-user=tu.mail@gmail.com
#SBATCH --output=./logs/%A_%a_%x_%j.out
#SBATCH --cpus-per-task=1
#SBATCH --array=1-2

# Ruta aos xenomas de referencia
# AQUÍ O QUE DAMOS É A RUTA ONDE ESTÁN TODAS AS POBOACIÓNS
REF_DIR="/mnt/lustre/scratch/nlsas/home/usc/gb/sdd/PRSs_IMPORTANTE/prscsx/ref_pops/"

# Prefixo do arquivo BIM (sen ".bim")
BIM_PREFIX="/ruta/bim_rs_allchr"

# Summary statistics formateados de todas as POP separados por coma
SUMSTATS="/ruta/sumstats/formateados/POP1.txt,/ruta/sumstats/formateados/POP2.txt"

# N do GWAS (tamaños de mostra, ordeados por POP)
N_GWAS="NUMERO_POP1,NUMERO_POP2"

# Cromosoma que estamos a analizar
CHR="1"

# Poboacións (EUR,AMR,AFR,SAS,EAS), mesma orde que SUMSTATS
POPS="POP1,POP2"

# Directorio de saída
OUT_DIR="/ruta/output"

# Nome do ficheiro de saída
OUT_NAME="nome_saida_prscsx"

# Opción META (True ou False)
META="True"

# Cargar entorno virtual
module load miniconda3/4.11.0
source activate prscsx

# Executar PRScsx
python /mnt/lustre/scratch/nlsas/home/usc/gb/sdd/PRSs_IMPORTANTE/prscs_prscsx/PRScsx//PRScsx.py \ 
  --ref_dir="${REF_DIR}" \
  --bim_prefix="${BIM_PREFIX}" \
  --sst_file="${SUMSTATS}" \
  --n_gwas="${N_GWAS}" \
  --chrom="${CHR}" \
  --pop="${POPS}" \
  --out_dir="${OUT_DIR}" \
  --out_name="${OUT_NAME}" \
  --meta="${META}"

STATUS=$?
if [ "$STATUS" -eq 0 ]; then
    echo "CHR${CHR} SUCCESS" >> "${OUT_DIR}/${OUT_NAME}_prscsx_run_status.log"
else
    echo "CHR${CHR} FAILED with exit code $STATUS" >> "${OUT_DIR}/${OUT_NAME}_prscs_run_status.log"
fi
```

Nota: Podemos paralelizar o script, substituíndo o `$CHR` por `CHR=${SLURM_ARRAY_TASK_ID}` e especificando o número de arrays:

`#SBATCH --output=./logs/%A_%a_%x_%j.out
#SBATCH --cpus-per-task=1
#SBATCH --array=1-22`

⚠POSIBLE ERRO: 

Que exista algún beta/SE con NA/INF, ou OR raras. Se os *sumstats* non son moi fiables (N pequena), revisar os mínimos e máximos destas 3 variables, e detectar valores moi extremos ou inviables como OR=0. 

## Paso 3. Aplicamos aos datos target

Unha vez obtemos os *score* (un arquivo por cromosoma), temos que facer un paso para recuperar os SNP ID (recordemos que PRS-cs traballa con *rsid*). 

Primeiro, xuntamos todos os cromosomas coa liña seguinte (debemos cambiar o `${out_name}`co que queramos):

```bash
$out_name="output_name"

cat ${out_name}_pst_eff_a1_b0.5_phiauto_chr* > scores_PRScs_${out_name}_allchr.txt
```

**NOTAS**: No caso de PRS-CSx, os que debemos concatenar son os “META”. PRS-cs non ten un parámetro de OUT_NAME, pero si de OUT_DIR no que podemos especificar o prefixo.

### 01. Recuperamos os SNP-ID

O PRS aplícase sobre os datos xenéticos da cohorte *target*, que ten os nomes dos SNPs no `.bim`como `chr:pos:ref:alt`. Se aplicáramos os *scores* tal como saen de PRS-cs, non coincidiría ningunha variante. Debemos, por tanto, substituír os *rsid* cos códigos de SNP tipo `chr:pos:ref:alt`  no arquivo `scores` . 

Ademais, recordemos que non todas as variantes tiñan un equivalente tipo *rsid*: debemos excluílas dos arquivos binarios de `plink`. Se ben a función `--score` só usa os SNPs do arquivo de `scores`, imos máis seguras extraendo os SNPs e facendo un único binario con todos os cromosomas (explicado máis adiante).

Anteriormente creamos un arquivo chamado `equivalencia.rs.snp...allchr.txt` , e será a partir deste do que recuperemos os SNP ID.

***Nota***: O output de PRS-cs e PRS-csx ten as columnas CHR RSID POS A1 (effect allele) A2 (non-effect allele) EFFECT_SIZE (posterior). 

Usamos o seguinte script de R:

```r
rsid_conv <- read.table("ruta/ao/equivalencia.rsid.snp.txt", 
                        header = TRUE)
names(rsid_conv) <- c("rsid", "SNP")

# Abrimos o score dos PRS (scores_PRScs_${out_name}_allchr.txt)
scores <- read.table("scores_PRScs_OUTNAME_allchr.txt", header = TRUE)
names(scores) <- c("CHR", "rsid", "pos", "EA", "NEA", "score")

# Merge de scores con identificadores SNP
scores_snp <- merge(rsid_conv, scores, by = "rsid")

# Exportar resultado final
write.table(scores_snp, "NOME_OUTPUT_SCORES.txt", quote = FALSE, row.names = FALSE, sep = "\t")
```

Despois, exportamos unha lista dos SNPs incluídos no PRS para extraelos dos binarios:

```bash
cut -f2 -d$'\t' "/mnt/lustre/scratch/nlsas/home/usc/gb/sdd/PRS_TESIS/prscsx/prscsx_only_eur/scores_snp_prscs_eur_noSC.txt" > extract_prscs_eur.txt
```

### 02. Extraemos os SNPs dos datos target e combinamos os chr.

Aquí, a estratexia seguida é a de seleccionar os SNPs do PRS, extraelos dos binarios *target*, e facer un merge deses binarios con todos os cromosomas. Outra estratexia sería facer primeiro o *merge* e despois a extracción dos SNPs (menos arquivos temporais e menos espacio en disco), pero é mellor como está abaixo porque ao facer o *merge* pode haber conflitos con SNPs multialélicos, etc. A lista de PRS-cs xa está curada e, salvo excepcións, non acarrea este tipo de problemas. 

**NOTA:** Aquí volvemos a asumir que se seguiron os nosos pasos da imputación. Do contrario, habería que modificar o `--bfile` co nome correspondente (aquí, por exemplo, temos o sufixo “_imputado”. Como borramos os arquivos temporais por defecto para que non ocupen espazo, é importante que exportemos os erros.

```bash
#!/bin/bash

module load plink

# 1. Configuración de rutas e nomes
TARGET_DATA_FOLDER="ruta/a/carpeta/con/binarios/imputados"
INPUT_DATA_NAME="nome_binarios_sin_chr"
SNPS_EXTRACT="archivo/snps/extraer/para/score"
SCORE_FILE="ruta/a/score/file"
OUTPUT_FOLDER="ruta/a/output"
OUTPUT_SCORE_NAME="nome_de_scores_aplicados"

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

```

⚠NOTA IMPORTANTE: 

Os scores que usa PRS-Cs ou PRS-CSx van usar só as variantes dispoñibles no `.bim` para a construción do PRS, que á súa vez partirán das de HapMap (este método de construción de PRS só usa SNPs de HapMap). Ao ser genome-wide que coincidan exactamente coas do summary non é tan importante).
