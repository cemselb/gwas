plot_highlighted_manhattan <- function(gwas_results, known_loci, window_size = 100000, pval_threshold = 1e-1, title = "GWAS Results") {
  
  required_cols <- c("CHR", "BP", "P", "SNP")
  if (!all(required_cols %in% colnames(gwas_results))) {
    stop("gwas_results must contain CHR, BP, P, and SNP columns.")
  }
  
  clrs <- c("#6e0026", "#78e372", "#8d1c90", "#98e793", "#1f49ba", "#f1a930", "#001e6c", 
            "#ead475", "#df8aff", "#547100", "#ee6bd4", "#008143", "#d52357", "#62a1ff", 
            "#ff7649", "#453c7b", "#934c00", "#ffa3e5", "#632400", "#ff7d73", "#d1323c", "#ed9e79")
  
  filtered_results <- gwas_results[gwas_results$P < pval_threshold, ]
  
  known_merged <- merge(known_loci, filtered_results, by = "SNP", all.x = TRUE)
  known_merged <- known_merged[!is.na(known_merged$CHR), ]
  

  find_regional_snps <- function(chr, bp) {
    filtered_results$SNP[
      filtered_results$CHR == chr & 
      filtered_results$BP >= (bp - window_size) & 
      filtered_results$BP <= (bp + window_size)
    ]
  }
  
  snps_to_highlight <- unique(unlist(
    Map(find_regional_snps, known_merged$CHR, known_merged$BP)
  ))
  
  gws_line <- -log10(5e-8)
  
  qqman::manhattan(filtered_results, 
                   chr = "CHR", bp = "BP", p = "P", snp = "SNP", 
                   col = clrs, 
                   chrlabs = NULL, 
                   suggestiveline = 5, 
                   genomewideline = gws_line, 
                   highlight = snps_to_highlight,
                   ylim = c(0, 15),
                   main = title)
}
