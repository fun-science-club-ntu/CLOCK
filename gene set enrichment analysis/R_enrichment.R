library(clusterProfiler)
library(enrichplot)
library(plyr)
library(ggrepel)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

selectedGeneID <- c("CLOCK")

gsym.fc <- read.table("~/downloads/correlation_CLOCK_gene.txt", header = T,sep="\t")
gsym.fc <- gsym.fc[,c(1,4)]
colnames(gsym.fc)<-c("SYMBOL","logFC")
gsym.fc<-na.omit(gsym.fc)

mycol <- c("darkgreen","chocolate4","blueviolet","#223D6C","#D20A13","#088247","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")


#transfer gene symbol to ENTREZ ID
gsym.id <- bitr(gsym.fc$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
#head(gsym.id)
#dim(gsym.id)

gsym.fc.id <- merge(gsym.fc, gsym.id, by="SYMBOL", all=F)
#head(gsym.fc.id)

#rank the spearman rho estimate
gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T),]
#head(gsym.fc.id.sorted)

#let gene ID and spearman rho estimate as input
id.fc <- gsym.fc.id.sorted$logFC
names(id.fc) <- gsym.fc.id.sorted$ENTREZID
#head(id.fc)


#Geneset enrichment analysis of KEGG 
kk <- gseKEGG(id.fc, organism = "hsa")

dim(kk)
head(kk)

#transfer ENTREZ ID to gene symbol
kk.gsym <- setReadable(kk, 'org.Hs.eg.db', #Species, Human
                       'ENTREZID')

#rank pathways according to enrichment score
sortkk <- kk.gsym[order(kk.gsym$enrichmentScore, decreasing = T),]
head(sortkk)
tail(sortkk)

#write the enriched pathways
write.csv(sortkk,"gsea_output_CLOCK.csv", quote = F, row.names = F)

#
geneSetID <- c("hsa04012","hsa04350","hsa04120")

# select the gene of interest
selectedGeneID <- c("CLOCK")



#save the enrichment score of each pathway 
for (i in geneSetID) {
  gseaplot(kk, i)
  myGeneList <- enrichplot:::gsInfo(kk, i)
  row.names(myGeneList) <- gsym.fc$gsym
  myGeneList$id <- gsym.fc$ENTREZID 
  write.csv(myGeneList, paste0("gsea_genelist_", i, "_group1.csv"))
}




x <- kk
geneList <- position <- NULL ## to satisfy codetool

#combine multiple pathways 
gsdata <- do.call(rbind, lapply(geneSetID, enrichplot:::gsInfo, object = x))
gsdata$gsym <- rep(gsym.fc.id.sorted$SYMBOL,3)
#gsdata$gsym <- gsym.fc.id.sorted$SYMBOL



# 画running score
p.res <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
  geom_line(aes_(y = ~runningScore, color= ~Description), size=1) +
  scale_color_manual(values = mycol) +
  
  #scale_x_continuous(expand=c(0,0)) + #
  geom_hline(yintercept = 0, lty = "longdash", lwd = 0.2) + #
  ylab("Enrichment\n Score") +
  
  theme_bw() +
  theme(panel.grid = element_blank()) + #
  
  theme(legend.position = "top", legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent")) +
  
  theme(axis.text.y=element_text(size = 12, face = "bold"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank(),
        plot.margin=margin(t=.2, r = .2, b=0, l=.2, unit="cm"))
#p.res


# plot the ranking
rel_heights <- c(1.5, .5, 1.5) 

i <- 0
for (term in unique(gsdata$Description)) {
  idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
  gsdata[idx, "ymin"] <- i
  gsdata[idx, "ymax"] <- i + 1
  i <- i + 1
}
#head(gsdata)

p2 <- ggplot(gsdata, aes_(x = ~x)) +
  geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
  xlab(NULL) + ylab(NULL) + 
  scale_color_manual(values = mycol) + #use my favourite color
  
  theme_bw() +
  theme(panel.grid = element_blank()) + 
  
  theme(legend.position = "none",
        plot.margin = margin(t=-.1, b=0,unit="cm"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line.x = element_blank()) +
  scale_y_continuous(expand=c(0,0))
#p2


# plot the rho estimate
df2 <- p.res$data
df2$y <- p.res$data$geneList[df2$x]
df2$gsym <- p.res$data$gsym[df2$x]
#head(df2)

# plot the gene of interest
selectgenes <- data.frame(gsym = selectedGeneID)
selectgenes <- merge(selectgenes, df2, by = "gsym")
selectgenes <- selectgenes[selectgenes$position == 1,]
head(selectgenes)

p.pos <- ggplot(selectgenes, aes(x, y, fill = Description, color = Description, label = gsym)) + 
  geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0), 
               color = "grey") +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = mycol, guide=FALSE) + 
  scale_color_manual(values = mycol, guide=FALSE) + 
  
  #scale_x_continuous(expand=c(0,0)) +
  geom_hline(yintercept = 0, lty = 2, lwd = 0.2) + 
  ylab("Ranked list\n metric") +
  xlab("Rank in ordered dataset") +
  
  theme_bw() +
  theme(axis.text.y=element_text(size = 12, face = "bold"),
        panel.grid = element_blank()) +
  
  # plot the gene of interest
  geom_text_repel(data = selectgenes, 
                  show.legend = FALSE, #no figure legends
                  direction = "x", 
                  ylim = c(2, NA), 
                  angle = 90, 
                  size = 2.5, box.padding = unit(0.35, "lines"), 
                  point.padding = unit(0.3, "lines")) +
  theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))
p.pos


# make those plots together
plotlist <- list(p.res, p2, p.pos)
n <- length(plotlist)
plotlist[[n]] <- plotlist[[n]] +
  theme(axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.text.x = element_text(size = 12, face = "bold"))

plot_grid(plotlist = plotlist, ncol = 1, align="v", rel_heights = rel_heights)

# save the plot
ggsave("GSEA_multi_pathways_CLOCK.pdf", width=6, height=5)

