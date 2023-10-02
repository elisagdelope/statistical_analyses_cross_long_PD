# Title: func_plot_compareDEGs.R
# Authorship: Elisa Gomez de Lope, Contact details: elisa.gomezdelope@uni.lu
# Info: This script contains functions to compare 2 DE analysis results based on fdr or fdr + logfoldchange

library(ggplot2)

### ----------------------------------------------
# function to compare 2 DE analysis results, returns significant values based on false discovery rate (fdr)
# usage: 
# compare_DEGs(r1 = RES[["STAR"]],
#              r2 = RES[["salmon"]],
#              r1.lab = "STAR",
#              r2.lab = "Salmon",
#              out.file = "compare_DEGs.txt",
#              out.plot = "compare_DEGs_plot.pdf")

compare_DEGs_fdr <- function(r1, r2, r1.lab="condition1", r2.lab="condition2", alpha=0.05, out.file="compare_DEGs.txt", out.plot="compare_DEGs_plot.pdf") {
  ix <- match(rownames(r1), rownames(r2))
  r2 <- r2[ix, ]
  flt <- !is.na(r1$padj) & !is.na(r2$padj)
  r1 <- r1[flt, ]
  r2 <- r2[flt, ]
  df <- data.frame(
    gene = r1$gene_id,
    r1 = -log10(r1$padj),
    r2 = -log10(r2$padj)
  )
  
  pdf(file = out.plot)
  p <- ggplot(data = df, aes(x = r1, y = r2)) +
    geom_point(alpha = 0.5) +
    theme_bw() +
    xlab(r1.lab) +
    ylab(r2.lab) +
    xlim(c(0, max(df[, c("r1", "r2")]))) +
    ylim(c(0, max(df[, c("r1", "r2")]))) +
    geom_hline(yintercept = -log10(alpha), color = "red", lty = 2) +
    geom_vline(xintercept = -log10(alpha), color = "red", lty = 2) +
    ggtitle("-log10(adjusted p-values)")
  print(p)
  dev.off()
  
  r1.deg <- r1$gene_id[r1$padj <= alpha]
  r2.deg <- r2$gene_id[r2$padj <= alpha]
  
  print(paste0("# ", r1.lab, ": ", length(r1.deg)))
  print(paste0("# ", r2.lab, ": ", length(r2.deg)))
  print(paste0("# in common: ", length(intersect(r1.deg, r2.deg))))
  print("In common:")
  print(intersect(r1.deg, r2.deg))
  print(paste0("Unique to ", r1.lab, ": "))
  print(setdiff(r1.deg, r2.deg))
  print(paste0("Unique to ", r2.lab, ": "))
  print(setdiff(r2.deg, r1.deg))
  write_lines(c(paste0("# ", r1.lab, ": ", length(r1.deg), "\t# ", r2.lab, ": ", length(r2.deg), "\t# in common: ", length(intersect(r1.deg, r2.deg)))), file = out.file)
  write_lines("In common:", file = out.file, append = TRUE)
  write_lines(str_c(intersect(r1.deg, r2.deg), collapse = ','), file = out.file, append = TRUE)
  write_lines(paste0("Unique to ", r1.lab, ": "), file = out.file, append = TRUE)
  write_lines(str_c(setdiff(r1.deg, r2.deg), collapse = ','), file = out.file, append = TRUE)
  write_lines(paste0("Unique to ", r2.lab, ": "), file = out.file, append = TRUE)
  write_lines(str_c(setdiff(r2.deg, r1.deg), collapse = ','), file = out.file, append = TRUE)
  
  return(intersect(r1.deg, r2.deg))
}



# function to compare 2 DE analysis results, returns significant values based on false discovery rate (fdr) AND logfoldchange
compare_DEGs_fdr_lfchange <- function(r1, r2, r1.lab="condition1", r2.lab="condition2", alpha=0.05, lfchange=1, out.file="compare_DEGs.txt", out.plot="compare_DEGs_plot.pdf") {
  ix <- match(rownames(r1), rownames(r2))
  r2 <- r2[ix, ]
  flt <- !is.na(r1$padj) & !is.na(r2$padj)
  r1 <- r1[flt, ]
  r2 <- r2[flt, ]
  r1$is_lfchange <- ((r1$log2FoldChange >= lfchange) | (r1$log2FoldChange <= -lfchange))
  r2$is_lfchange <- ((r2$log2FoldChange >= lfchange) | (r2$log2FoldChange <= -lfchange))
  
  df <- data.frame(
    gene = r1$gene_id,
    r1 = -log10(r1$padj),
    r2 = -log10(r2$padj),
    lfchange = ((r1$is_lfchange == TRUE) | (r2$is_lfchange == TRUE))
  )
  
  pdf(file = out.plot)
  p <- ggplot(data = df, aes(x = r1, y = r2, color = lfchange)) +
    geom_point(alpha = 0.6) +
    theme_bw() +
    xlab(r1.lab) +
    ylab(r2.lab) +
    xlim(c(0, max(df[, c("r1", "r2")]))) +
    ylim(c(0, max(df[, c("r1", "r2")]))) +
    geom_hline(yintercept = -log10(alpha), color = "red", lty = 2) +
    geom_vline(xintercept = -log10(alpha), color = "red", lty = 2) +
    ggtitle("-log10(adjusted p-values)") +
    scale_color_manual(name = "log2 fold change >=1 | <= -1", values = c("FALSE" = "black", "TRUE" = "blue")) +
    theme(legend.position = c(0.7, 0.9))
  print(p)
  dev.off()
  
  r1.deg <- r1$gene_id[(r1$padj <= alpha) & ((r1$log2FoldChange >= lfchange) | (r1$log2FoldChange <= -lfchange))]
  r2.deg <- r2$gene_id[(r2$padj <= alpha) & ((r1$log2FoldChange >= lfchange) | (r1$log2FoldChange <= -lfchange))]
  
  print(paste0("# ", r1.lab, ": ", length(r1.deg)))
  print(paste0("# ", r2.lab, ": ", length(r2.deg)))
  print(paste0("# in common: ", length(intersect(r1.deg, r2.deg))))
  print("In common:")
  print(intersect(r1.deg, r2.deg))
  print(paste0("Unique to ", r1.lab, ": "))
  print(setdiff(r1.deg, r2.deg))
  print(paste0("Unique to ", r2.lab, ": "))
  print(setdiff(r2.deg, r1.deg))
  write_lines(c(paste0("# ", r1.lab, ": ", length(r1.deg), "\t# ", r2.lab, ": ", length(r2.deg), "\t# in common: ", length(intersect(r1.deg, r2.deg)))), file = out.file)
  write_lines("In common:", file = out.file, append = TRUE)
  write_lines(str_c(intersect(r1.deg, r2.deg), collapse = ','), file = out.file, append = TRUE)
  write_lines(paste0("Unique to ", r1.lab, ": "), file = out.file, append = TRUE)
  write_lines(str_c(setdiff(r1.deg, r2.deg), collapse = ','), file = out.file, append = TRUE)
  write_lines(paste0("Unique to ", r2.lab, ": "), file = out.file, append = TRUE)
  write_lines(str_c(setdiff(r2.deg, r1.deg), collapse = ','), file = out.file, append = TRUE)
  
  return(intersect(r1.deg, r2.deg))
}


