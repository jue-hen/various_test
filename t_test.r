###多个不同分组，多个不同连续变量
###data表示分析的数据框， group表示不同的类别变量， var表示要分析的连续变量
t_text <- function(data, group, var){
    t_res <- as.data.frame(matrix(ncol = 3, nrow = length(group)*length(var)))
    colnames(t_res) <- c("group","var","p_val")
    library(Hmisc)
    num=1
    for (i in c(1:length(group))) {
      for (j in c(1:length(var))) {
        #print(j)
        tmp_dat <- data[,c(group[i],var[j])]
        tmp_dat <- na.omit(tmp_dat)
        ###tmp_dat <- tmp_dat %>% filter(!is.na(tar_gen[i])) %>% filter(!is.na(tm[j]))
        colnames(tmp_dat) <- c("group","value")
        tmp_dat$group <- factor(tmp_dat$group) 
          ## 判断是否有多个分组
        if (length(levels(tmp_dat$group))) == 1){
          next
        }
        re<-t.test(value ~ group, data = tmp_dat, exact=F)
        
        t_res[num,1] <- group[i]
        t_res[num,2] <- var[j]
        t_res[num,3] <-re$p.value
        
        num = num+1
      }
    }
    t_res <- na.omit(t_res)
    t_res$fdr <- p.adjust(t_res$p_val, method = "fdr")
    return(t_res)
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
  wt_mean <- apply(wt, 2, mean)
  wt_sd<- apply(wt, 2, sd)
              
  mut <- dat[dat$genotype == "MUT",var]
  mut <- apply(mut, 2, function(x)(impute(x,median)))
  mut_mean <- apply(mut, 2, mean)
  mut_sd <- apply(mut, 2, sd)
  
  res <- as.data.frame(matrix(nrow = 8, ncol = 3,0))
  rownames(res) <- var
  colnames(res) <- c("WT","MUT","p")
  for (i in c(1:length(var))){

    wt_tmp <- paste(round(wt_mean,2),round(wt_sd,2), sep = "±")
    
    mut_tmp <- paste(round(mut_mean,2),round(mut_sd,2), sep = "±")
    
    res[i,1] <- wt_tmp
    res[i,2] <- mut_tmp
    res[i,3] <- t.test(wt[,i], mut[,i])$p.value
  }
  return(res)
}
