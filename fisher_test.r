###fisher 检验
library(tidyverse)

fish_test_calculate <- function(data, var1, var2){
    fish_res <- as.data.frame(matrix(nrow = length(var1)*length(var2), ncol = 3))
    colnames(fish_res) <- c("var1","var2","p","OR")
      
    num=1
    for (i in c(1:length(var1))) {
      for (j in c(1:length(var2))){
        ##fisher检验
        tmp_dat <- data.frame(v1=data[,var1[i]], v2=data[,var2[j]])
        tmp_dat <- tmp_dat %>% filter(!is.na(tmp_dat[,1])) %>% filter(!is.na(tmp_dat[,2]))
        
        tmp_dat <- as.data.frame(apply(tmp_dat, 2, as.character)) %>% lapply(tmp_dat,factor)
        ## 判断是否有多个分组
        if (length(levels(tmp_dat[,1])) == 1 | length(levels(tmp_dat[,2]))) == 1){
          next
        }

        #####进行检验，可以修改为卡方检验
        f_test <- fisher.test(table(tmp_dat[,1], tmp_dat[,2]))
        ##提取p值
        p_val <- f_test$p.value
        # formatted_p_val <- ifelse(p_val < 0.001, "<0.001", round(p_val, 3))
        ##提取OR值
        OR_val <- as.numeric(f_test$estimate)
    
        fish_res[num,1] <- var1[i]
        fish_res[num,2] <- var2[j]
        fish_res[num,3] <- p_val
        fish_res[num,4] <- OR_val
      }
      num = num + 1
    }
    
    fish_res$fdr <- p.adjust(fish_res$p, method = "fdr"))
    fish_res <- fish_res[,c(1,2,3,5,4)]
    return(fish_res)
}
