library(data.table)
library(Gviz)
library(RColorBrewer)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
txdb_mm9 <- TxDb.Mmusculus.UCSC.mm9.knownGene
grt <- GeneRegionTrack(txdb_mm9, genome="mm9",showId=TRUE, geneSymbol=TRUE, name="UCSC")
library(org.Mm.eg.db)
z <- mapIds(org.Mm.eg.db, gene(grt), "SYMBOL", "ENTREZID", multiVals = "first")
zz <- sapply(z, is.null)
z[zz] <- gene(grt)[zz]
gr <- ranges(grt)
mcols(gr)$symbol <- z
grt@range <- gr

bwInfo<-read.table("bwfilelist",header=F,row.names=1,as.is=T)
head(bwInfo)


startpoint<-100000
endpoint<-200000
chr<-"chr5"


tracklist<-list()
itrack <- IdeogramTrack(genome = "mm9", chromosome = chr,outline=T)
tracklist[["itrack"]]<-itrack
scalebar <- GenomeAxisTrack(scale=0.25,col="black",fontcolor="black",name="Scale",labelPos="above",showTitle=TRUE)
tracklist[["scalebar"]]<-scalebar
axisTrack <- GenomeAxisTrack(labelPos="above",col="black",fontcolor="black",name=paste(chr,":",sep=""),exponent=0,showTitle=TRUE)
tracklist[["axisTrack"]]<-axisTrack
colpal<-rep(brewer.pal(12,"Dark2"),20)
coldf<-data.frame(col=colpal[1:nrow(bwInfo)],row.names = rownames(bwInfo),stringsAsFactors = F)
for(index in rownames(bwInfo)){
  bgFile<-file.path("please input the bw file path",paste(index,".bw",sep=""))
  tracklist[[index]]<-DataTrack(range = bgFile,genome="mm9",type="histogram",
                                name=chartr("_","\n",bwInfo[index,]),
                                col.histogram=coldf[index,])
}
tracklist[["grt"]]<-grt
plotTracks(tracklist, from = startpoint, to = endpoint,
           chromosome=chr,background.panel = "white", background.title = "white",
           col.title="black",col.axis="black",
           rot.title=0,cex.title=0.9,margin=38,title.width=1.5,
           collapseTranscripts = "longest",
           ylim = c(0, 30))