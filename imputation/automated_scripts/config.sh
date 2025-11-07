

#########################################################
#-------------- COMO USAR -------------------------------
#########################################################

# Modificamos este script segundo indicamos abaixo
# Despois executamos Script_pre_imputación.sh 
# Os arquivos que subimos a imputar son os chr${i}_flipped_preimputacion.vcf.gz
# E unha vez os teñamos imputados, Script_post_imputación.sh
# Os que usaremos para gwas, PRS, etc. son os ${PREFIX_IMPUTADOS}_chr${i}_imputado



#----------------------------------------------------------------------------
#-----# MODIFICAR ANTES E DESPOIS DA IMPUTACIÓN SEGUNDO SE INDICA -----------
#----------------------------------------------------------------------------


#---------- A CAMBIAR PRE-IMPUT

DIR_XENOT="/ruta/directorio/xenotipados"
DIR_PRE_IMPUTACION="/ruta/onde/queramos/gardar/arquivos/pre-imputación"
# podemos usar a memsa ruta que en DIR_XENOTIPADOS
# simplemente o separo por se queremos dúas carpetas diferentes

NOME_XENOT="nombre_arquivo_binario_dos_xenotipados"
GENOME_BUILD="hg19" 
# ou hg38

#--------------- NON CAMBIAR -------------------------------------------------
FASTA_hg19="/mnt/netapp1/Store_chumxrcg/REFERENCE_GENOMES/hg19.fa"
FASTA_hg38="/mnt/netapp1/Store_chumxrcg/REFERENCE_GENOMES/hg38.fa"
RENAME_FILE="/mnt/netapp1/Store_chumxrcg/REFERENCE_GENOMES/RENAME_hg19_back.txt" 
#-----------------------------------------------------------------------------


#---------- A CAMBIAR POST-IMPUT

PASSWORD="contrasinal_imputación"
LINK_DESCARGA="link_que_da_topmed_para_descarga"
DIR_IMPUTACION="/ruta/directorio/para/gardar/imputados"
PREFIX_IMPUTADOS="nome_saída_imputados"
MAF="0.01"
