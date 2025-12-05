######## Lánzase con 
# source("anotar_bim_con_rsid_topmed.R")
# anotar_bim_con_rsid_topmed(chr=1:23,bim_path_prefix = "/mnt/netapp1/Store_chumxrcg/EB/EB25/02_imputacion_eur0725/imputados/EB25_EUR_", 
# 										info_path = "/mnt/netapp1/Store_chumxrcg/EB/EB25/02_imputacion_eur0725/imputados/",
# 										out_path = "/mnt/netapp1/Store_chumxrcg/EB/EB25/02_imputacion_eur0725/creacion_bim_rsid_para_prscs/", 
# 										out_prefix = "EBeur")


anotar_bim_con_rsid_topmed <- function(chr = 1:23,
                                 bim_path_prefix = "/ruta/a/binarios/imputados/prefixo",
                                 info_path = "/ruta/a/info/",
																 out_path = "/ruta/a/carpeta/prs/",
                                 out_prefix = "nome_output") {
	 if (!dir.exists(out_path)) {
    dir.create(out_path, recursive = TRUE)
	}
  
  for (i in chr) {

    # 1. Ler o  BIM do cromosoma i
    ficheiro_bim <- paste0(bim_path_prefix, "chr", i, "_imputado.bim")
    bim <- read.table(ficheiro_bim)
    names(bim) <- c("CHR", "snp", "cm", "POS", "A1", "A2")
    
   # 2. Lemos o INFO
    if(i == 23){
    ficheiro_info <- paste0(info_path,"/chrX.info.gz")
    } else {
    ficheiro_info <- paste0(info_path, "/chr", i, ".info.gz")}
    
    info <- read.table(ficheiro_info)
    names(info) <- c("CHR", "POS", "RSID", "REF", "ALT", "QUAL", "FILTER", "INFO")

    # 3. Crear SNPID no formato chr:pos:ref:alt
    info$SNP <- paste0(info$CHR, ":", info$POS, ":", info$REF, ":", info$ALT)
    
    # 4. Facer merge co BIM a través do SNPID
    bim_rs <- merge(bim, info[c("RSID", "SNP")], by.x = "snp", by.y = "SNP", sort = FALSE)
    
    # 5. Filtrar SNPs rsid válido
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

