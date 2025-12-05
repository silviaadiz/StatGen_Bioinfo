##### Lanzado dende bash con:
#   Rscript --vanilla "/mnt/lustre/scratch/nlsas/home/usc/gb/sdd/EB25/PRSs_novembro25_cohortes_combi/PRS_psy/Apply_to_EB/2_Convert_rsid_to_SNP_opt.R" \
#   --equivalence_file "/mnt/lustre/scratch/nlsas/home/usc/gb/sdd/EB25/PRS_Docherty/00_archivos_input_PRS/BIM_rs_PRScs/equivalencia.rs.snp.EBeur_senX.txt" \
#   --score "${SCOREFILE_OUTPUT}" --pheno_name ${pheno} --target_name ${TARGET_NAME}



library(optparse)
option_list=list(
  make_option(c("--equivalence_file"),type="character"),
  make_option(c("--score"),type="character"),
  make_option(c("--pheno_name"),type="character"),
  make_option(c("--target_name"),type="character"))
  
opt=parse_args(OptionParser(option_list=option_list))
 


rsid_conv <- read.table(opt$equivalence_file, header = TRUE)
names(rsid_conv) <- c("rsid", "SNP")

# Abrimos o score dos PRS (scores_PRScs_${out_name}_allchr.txt)
scores <- read.table(opt$score, header = TRUE)
# É o PRS aplicado a europeos, así que poñemos o bim de EUR!!!!!


names(scores) <- c("CHR", "rsid", "pos", "EA", "NEA", "score")

# Merge de scores con identificadores SNP
scores_snp <- merge(rsid_conv, scores, by = "rsid")

# Exportar resultado final
outp_name <- paste0("SNP.RSID.Scores_prscs_",opt$pheno_name,"_to_",opt$target_name,".txt")
write.table(scores_snp, outp_name, quote = FALSE, row.names = FALSE, sep = "\t")
