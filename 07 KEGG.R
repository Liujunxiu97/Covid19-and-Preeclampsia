#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")
#install.packages("circlize")
#install.packages("RColorBrewer")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")
#BiocManager::install("ComplexHeatmap")
#Import packages.
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)

pvalueFilter=0.05      #p-value filtering criteria.
qvalueFilter=0.05      #Filtering criteria for adjusted p-values.

#Define colors.
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}
	
setwd("XXXXX")          #Set the working directory.
rt=read.table("interGenes.txt", header=T, sep="\t", check.names=F)     #Read the input file.

#Extract the names of differentially expressed genes and convert gene names to gene IDs.
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #Remove genes with gene IDs that are 'NA'.
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#Perform KEGG enrichment analysis.
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#Save the significantly enriched results.
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#Define the number of displayed pathways.
showNum=15
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

#Bar chart.
pdf(file="barplot1.pdf", width=9, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, color=colorSel)
dev.off()

#Bubble plot.
pdf(file="bubble1.pdf", width = 9, height = 7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, color=colorSel)
dev.off()

###########Plot a KEGG circle plot.###########
Pathway.col=c("#90EE90", "#E7B800", "#00AFBB")
showNum=15
data=KEGG[order(KEGG$p.adjust),]
if(nrow(KEGG)>showNum){
	data=data[1:showNum,]
}
data$Pathway="KEGG"
main.col = Pathway.col[as.numeric(as.factor(data$Pathway))]

#Organize the data for the circle plot.
BgGene = as.numeric(sapply(strsplit(data$BgRatio,"/"),'[',1))
Gene = as.numeric(sapply(strsplit(data$GeneRatio,'/'),'[',1))
ratio = Gene/BgGene
logpvalue = -log(data$pvalue,10)
logpvalue.col = brewer.pal(n = 8, name = "Reds")
f = colorRamp2(breaks = c(0,2,4,6,8,10,15,20), colors = logpvalue.col)
BgGene.col = f(logpvalue)
df = data.frame(KEGG=data$ID,start=1,end=max(BgGene))
rownames(df) = df$KEGG
bed2 = data.frame(KEGG=data$ID,start=1,end=BgGene,BgGene=BgGene,BgGene.col=BgGene.col)
bed3 = data.frame(KEGG=data$ID,start=1,end=Gene,BgGene=Gene)
bed4 = data.frame(KEGG=data$ID,start=1,end=max(BgGene),ratio=ratio,col=main.col)
bed4$ratio = bed4$ratio/max(bed4$ratio)*9.5

#Plot the main body of the circle plot.
pdf(file="KEGG.circlize.pdf",width=10,height=10)
par(omi=c(0.1,0.1,0.1,1.5))
circos.par(track.margin=c(0.01,0.01))
circos.genomicInitialize(df,plotType="none")
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA,bg.col = main.col)

for(si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si,track.index = 1,
              major.at=seq(0,max(BgGene),by=100),labels.facing = "clockwise")
}
f = colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
circos.genomicTrack(bed2, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[,2], 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], adj=0,cex=0.8,...)
                    })
circos.genomicTrack(bed3, ylim = c(0, 1),track.height = 0.1,bg.border="white",
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], cex=0.9,adj=0,...)
                    })
circos.genomicTrack(bed4, ylim = c(0, 10),track.height = 0.35,bg.border="white",bg.col="grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim = get.cell.meta.data("cell.xlim")
                      cell.ylim = get.cell.meta.data("cell.ylim")
                      for(j in 1:9) {
                        y = cell.ylim[1] + (cell.ylim[2]-cell.ylim[1])/10*j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], col = value[,2], 
                                         border = NA, ...)
                      #circos.genomicText(region, value, y = 0.3, labels = value[,1], ...)
                    })
circos.clear()
#Plot the legend in the middle of the circle plot.
middle.legend = Legend(
  labels = c('Number of Genes','Number of Select','Rich Factor(0-1)'),
  type="points",pch=c(15,15,17),legend_gp = gpar(col=c('pink','#BA55D3',Pathway.col[1])),
  title="",nrow=3,size= unit(3, "mm")
)
circle_size = unit(1, "snpc")
draw(middle.legend,x=circle_size*0.42)
#Plot the legend for KEGG categories.
main.legend = Legend(
  labels = c("KEGG"),  type="points",pch=15,
  legend_gp = gpar(col=Pathway.col), title_position = "topcenter",
  title = "Pathway", nrow = 3,size = unit(3, "mm"),grid_height = unit(5, "mm"),
  grid_width = unit(5, "mm")
)
#Plot the legend for p-values.
logp.legend = Legend(
  labels=c('(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,15]','(15,20]','>=20'),
  type="points",pch=16,legend_gp=gpar(col=logpvalue.col),title="-log10(Pvalue)",
  title_position = "topcenter",grid_height = unit(5, "mm"),grid_width = unit(5, "mm"),
  size = unit(3, "mm")
)
lgd = packLegend(main.legend,logp.legend)
circle_size = unit(1, "snpc")
print(circle_size)
draw(lgd, x = circle_size*0.85, y=circle_size*0.55,just = "left")
dev.off()


