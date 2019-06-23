tcga_expr<-read.table("tcga_gene_tpm_filtercancers.txt",header=T,sep="\t",row.names=1)
#tcga_RSEM_gene_tpm_filtercancers.txt, this file is huge. Download link is provided in our paper.


library(reshape2)
library(plyr)
library(ggplot2)
symbol <- "ENSG00000069667.15" #RORA, you can replace this gene to your interested one

# target gene
target.exps <- tcga_expr[symbol,]
# other gene
other.expr <- tcga_expr[-which(rownames(tcga_expr)==symbol),]

sgcor <- as.data.frame(cor(t(other.expr), t(target.exps))) 
colnames(sgcor) <- "r_pearson"
sgcor$pval_pearson <- apply(other.expr, 1, function(x) (cor.test(x, t(target.exps),na.omit=T)$p.value))
sgcor$r_pearson <- apply(other.expr, 1, function(x) (cor.test(x, t(target.exps),na.omit=T)$estimate))
sgcor_spearman <- as.data.frame(cor(t(other.expr), t(target.exps), method = "spearman"))
colnames(sgcor_spearman) <- "r_spearman"
sgcor_spearman$pval_spearman <- apply(other.expr, 1, function(x)(cor.test(x, t(target.exps), method = "spearman")$p.value))
sgcor_spearman$r_spearman <- apply(other.expr, 1, function(x)(cor.test(x, t(target.exps), method = "spearman")$estimate))
cors <- cbind(sgcor, sgcor_spearman)
cors$gene <- rownames(other.expr)
head(cors)

write.table(cors,"correlation_RORA.txt",quote=F,sep="\t")
