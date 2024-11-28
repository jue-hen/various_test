###exp为数据框， target为想要分析的变量， var为想要与目标变量进行相关性分析的变量， method表示相关分析方法，默认为皮尔逊

coefficient_test <- function(exp, target, var, method="pearson"){
    res <- as.data.frame(matrix(nrow = length(var), ncol = 3, 0))
    colnames(res) <- c("var","coefficient","p_value")
    
    for (i in c(1:length(var))) {
      gene_id <- var[i]
      tmp <- cor.test(exp[,target], exp[,gene_id], method = ) ###"pearson" "spearman"
      res[i,1] <- gene_id
      res[i,2] <- tmp$estimate
      res[i,3] <- tmp$p.value
    }
    res$fdr <- p.adjust(res$p_value, method = "fdr")
    return(res)
}
