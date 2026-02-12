expr_file <- "D:/out.csv"
deg_file <- "D:/deg.csv"
top10_output <- "D:/each_deg_top10_corr_genes.csv"
cor_matrix_output <- "D:/all_genes_deg_cor_matrix.csv" 


expr <- read.csv(expr_file, stringsAsFactors = FALSE, check.names = FALSE)
rownames(expr) <- expr$gene_id
expr_values <- expr[, colnames(expr) != "gene_id"]
for (i in 1:ncol(expr_values)) {
  expr_values[, i] <- as.numeric(as.character(expr_values[, i]))
}
expr_values <- expr_values[!apply(is.na(expr_values), 1, all), ]

deg <- read.csv(deg_file, stringsAsFactors = FALSE)
deg_genes <- deg[, 1]
matched_deg <- deg_genes[deg_genes %in% rownames(expr_values)]
if (length(matched_deg) == 0) stop("无匹配的DEG基因")
message(paste("✅ 匹配到", length(matched_deg), "个DEG基因"))


expr_t <- t(expr_values) 

full_cor_matrix <- cor(expr_t[, matched_deg], expr_t, method = "spearman")

full_cor_matrix <- t(full_cor_matrix)
cor_matrix_df <- data.frame(
  gene_id = rownames(full_cor_matrix),
  full_cor_matrix,
  stringsAsFactors = FALSE
)

all_top10 <- data.frame()
for (deg_gene in matched_deg) 
  cor_vec <- full_cor_matrix[, deg_gene]

  cor_vec_filtered <- cor_vec[!names(cor_vec) %in% deg_gene & !is.na(cor_vec)]
  cor_sorted <- sort(abs(cor_vec_filtered), decreasing = TRUE)

  top10 <- cor_sorted[1:min(10, length(cor_sorted))]
  deg_result <- data.frame(
  
    source_deg = deg_gene,

    correlated_gene = names(top10),

    spearman_cor_abs = as.numeric(top10),

    stringsAsFactors = FALSE
  )
  all_top10 <- rbind(all_top10, deg_result)
}


write.csv(cor_matrix_df, cor_matrix_output, row.names = FALSE)

write.csv(all_top10, top10_output, row.names = FALSE)

message(paste("完整相关性矩阵已输出至：", cor_matrix_output))
message(paste("每个DEG的Top10相关基因已输出至：", top10_output))
message(paste("相关性矩阵维度：", nrow(cor_matrix_df), "个基因 ×", length(matched_deg), "个DEG"))
message(paste("Top10列表共", nrow(all_top10), "条结果"))