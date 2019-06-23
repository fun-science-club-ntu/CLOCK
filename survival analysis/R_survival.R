library(survival)
library(survminer)

#format
#sample time  status  gene_1_expression gene_1_expression gene_x_expression
#sample_1 300 1 1.5 2.5 3.5
#sample_2 600 0 3.5 4.5 5.5
#sample_x 900 1 5.5 6.5 7.5


#define each gene cutoff
res.cut <- surv_cutpoint(svdata, time = "time", 
                         event = "status", 
                         variables = names(svdata)[3:ncol(svdata)], 
                         minprop = 0.05) 

res.cat <- surv_categorize(res.cut)
my.surv <- Surv(res.cat$futime, res.cat$fustat)
pl<-list()


#plot survival
for (i in colnames(res.cat)[3:ncol(svdata)]) {
  group <- res.cat[,i] 
  survival_dat <- data.frame(group = group)
  fit <- survfit(my.surv ~ group)
  group <- factor(group, levels = c("low", "high"))
  data.survdiff <- survdiff(my.surv ~ group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  svsort <- svdata[order(svdata[,i]),]
  pl[[i]]<-ggsurvplot(fit, data = survival_dat ,
                      conf.int = F, 
                      censor = F, 
                      palette = c("#B85A1F","#248868"),
                      break.x.by=300,
                      legend.title = i,
                      font.legend = 11,
                      legend.labs=c(paste0(">",round(svsort[fit$n[2],i],2),"(",fit$n[1],")"),
                                    paste0("<",round(svsort[fit$n[2],i],2),"(",fit$n[2],")")),
                      pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                 paste("p = ",round(p.val,3), sep = "")),
                                   HR, CI, sep = "\n"))
    ggsave(paste0(i,".pdf"),width = 6,height = 6)
}