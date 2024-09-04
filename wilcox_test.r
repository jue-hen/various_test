###基因型分型             
tar_gen <- colnames(analysis_dat)[c(6,7,9,10,11,20,21,22,24,25:220)]
tm <- c(dns_col,pct_col)
library(tidyverse)
analysis_dat <- analysis_dat %>%
  mutate_at(vars(dns_col), as.numeric)%>%
  mutate_at(vars(pct_col), as.numeric)

###多个不同分组，多个不同连续变量
###data表示分析的数据框， group表示不同的类别变量， var表示要分析的连续变量
wilcox_text <- function(data, group, var){
    wilcox_res <- as.data.frame(matrix(ncol = 3, nrow = length(group)*length(var)))
    colnames(wilcox_res) <- c("group","var","p_val")
    library(Hmisc)
    num=1
    for (i in c(1:length(group))) {
      for (j in c(1:length(var))) {
        #print(j)
        tmp_dat <- analysis_dat[,c(group[i],var[j])]
        tmp_dat <- na.omit(tmp_dat)
        #tmp_dat <- tmp_dat %>% filter(!is.na(tar_gen[i])) %>% filter(!is.na(tm[j]))
        colnames(tmp_dat) <- c("group","value")
        tmp_dat$group <- factor(tmp_dat$group) 
          ## 判断是否有多个分组
        if (length(levels(tmp_dat$group))) == 1){
          next
        }
        re<-wilcox.test(value ~ group, data = tmp_dat, exact=F)
        
        wilcox_res[num,1] <- group[i]
        wilcox_res[num,2] <- var[j]
        wilcox_res[num,3] <-re$p.value
        
        num = num+1
      }
    }
    wilcox_res <- na.omit(wilcox_res)
    wilcox_res$fdr <- p.adjust(wilcox_res$p_val, method = "fdr")
    return(wilcox_res)
}


###单个分组，多个不同的连续变量
###dat表示需要分析的数据框，数据框中需要含有genotype该列表示是否为突变个体， var表示需要分析的连续变量
library(Hmisc)
extract_dat <- function(dat, var){
  dat$genotype <- ifelse(dat$genotype == "0", "WT", "MUT")
  #table(dat$genotype)
  dat$genotype <- factor(dat$genotype, levels = c("WT","MUT"))
  
  wt <- dat[dat$genotype == "WT",var]
  wt <- apply(wt, 2, function(x)(impute(x,median)))
  wt_res <- apply(wt, 2, quantile)
  mut <- dat[dat$genotype == "MUT",var]
  mut <- apply(mut, 2, function(x)(impute(x,median)))
  mut_res <- apply(mut, 2, quantile)
  
  res <- as.data.frame(matrix(nrow = 8, ncol = 3,0))
  rownames(res) <- var
  colnames(res) <- c("WT","MUT","p")
  for (i in c(1:length(var))){
    wt_tmp1 <- wt_res[2,i]
    wt_tmp2 <- wt_res[4,i]
    wt_tmp <- paste(round(wt_tmp1,2),round(wt_tmp2,2), sep = "~")
    
    mut_tmp1 <- mut_res[2,i]
    mut_tmp2 <- mut_res[4,i]
    mut_tmp <- paste(round(mut_tmp1,2),round(mut_tmp2,2), sep = "~")
    
    res[i,1] <- wt_tmp
    res[i,2] <- mut_tmp
    res[i,3] <- wilcox.test(wt[,i], mut[,i])$p.value
  }
  return(res)
}
