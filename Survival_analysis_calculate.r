###生存分析提取p值，OR值等
###extract p 值
library("survival")
library("survminer")
    

Survival_analysis_calculate <- function(data, var){
    dat <- as.data.frame(matrix(nrow = length(var), ncol = 7))
    colnames(dat) <- c("var","pfs_p","pfs_95%CI","pfs_HR","os_p","os_95%CI","os_HR")
    
    for (j in c(1:length(genes))){
      i=genes[j]
      dat[j,1] <- i
      
      surv_object_pfs <- Surv(time=data$pfs_month, event = data$pfs_status)
      ###cox 比例风险模型
      cox_model_pfs <- coxph(surv_object_pfs ~ data[,i], data = data)
      ##提取HR和置信区间
      hr_pfs <- summary(cox_model_pfs)$coefficients[1,2]
      ci_lower_pfs <- summary(cox_model_pfs)$conf.int[1,3]
      ci_upper_pfs <- summary(cox_model_pfs)$conf.int[1,4]
      p_value_pfs <- summary(cox_model_pfs)$coefficients[1,5]
      dat[j,2] <- p_value_pfs
      dat[j,3] <- paste0(round(ci_lower_pfs,2),"-",round(ci_upper_pfs,2))
      dat[j,4] <-  round(hr_pfs,2)
    
      ###os
      surv_object_os <- Surv(time=data$os_month, event = data$os_status)
      ###cox 比例风险模型
      cox_model_os <- coxph(surv_object_os ~ data[,i], data = data)
      ##提取HR和置信区间
      hr_os <- summary(cox_model_os)$coefficients[1,2]
      ci_lower_os <- summary(cox_model_os)$conf.int[1,3]
      ci_upper_os <- summary(cox_model_os)$conf.int[1,4]
      p_value_os <- summary(cox_model_os)$coefficients[1,5]
      dat[j,5] <- p_value_os
      dat[j,6] <- paste0(round(ci_lower_os,2),"-",round(ci_upper_os,2))
      dat[j,7] <-  round(hr_os,2)
    }
    dat$pfs_fdr <- p.adjust(dat$pfs_p, method = "fdr")
    dat$os_fdr <- p.adjust(dat$os_p, method = "fdr")
    dat <- dat[,c(1,2,8,3,4,5,9,6,7)]
    return(dat)
}

