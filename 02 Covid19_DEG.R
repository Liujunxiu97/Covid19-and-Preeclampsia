rm(list = ls())
library("DESeq2")
library("limma")
library('futile.logger')
library( 'pheatmap' )
library( 'ggplot2' )
#install.packages('FactoMineR')
#install.packages('factoextra')
library('FactoMineR')
library(tinyarray)
library(ggpubr)
library(ggthemes)
#install.packages("flog.appender")
setwd("xxx")

padj = 0.05
logFoldChange=0.5

library(data.table)

expr<- read.table("mRNA.symbol.txt",sep = "\t" , header = T,
                  stringsAsFactors = FALSE ,
                  check.names = FALSE)

expr = avereps(expr[,-1],ID = expr$ENSEMBLID) 
expr = expr[rowMeans(expr)>1,]
#Expressing data organization.
#expr <- (2^expr - 1)
expr[1:5,1:5]
expr <- round(expr,0)
expr[1:5,1:5]
data=expr

#Read all files in the directory that end with "s1.txt".
sampleName1=c()
files=dir()
files=grep("s1.txt$", files, value=T)
for(file in files){
  rt=read.table(file, header=F, sep="\t", check.names=F)      #Read the input file.
  geneNames=as.vector(rt[,1])      #Extract gene names.
  uniqGene=unique(geneNames)       #Get unique genes.
  sampleName1=c(sampleName1, uniqGene)
}

#Read all files in the directory that end with "s2.txt".
sampleName2=c()
files=dir()
files=grep("s2.txt$", files, value=T)
for(file in files){
  rt=read.table(file, header=F, sep="\t", check.names=F)      #Read the input file.
  geneNames=as.vector(rt[,1])      #Extract gene names.
  uniqGene=unique(geneNames)       #Get unique genes.
  sampleName2=c(sampleName2, uniqGene)
}
#Extract data for the experimental group and the control group.
conData=data[,sampleName1]
treatData=data[,sampleName2]
data=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)
group_list <- c(rep('con',conNum),rep('treat',treatNum))

#Perform differential expression analysis using the DESeq2 package.
condition = factor(group_list)
coldata <- data.frame(row.names = colnames(data), condition)
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = coldata,
                              design = ~condition)
dds$condition<- relevel(dds$condition, ref = "con") # Specify which group to be designated as the control group.

#Differential expression matrix.
dds <- DESeq(dds)  
all <- as.data.frame(results(dds))
all = na.omit(all)
allout=all[,colnames(all)]
allout=rbind(ID=colnames(all),all)
write.table(allout,file="1.mRNA_all.txt",sep="\t",quote=F,col.names=F)

#Output corrected expression values.
Type=c(rep("con",conNum),rep("treat",treatNum))
outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData, file="normalize.txt", sep="\t", quote=F, col.names=F)


#Write protein and lncRNA table.
diffSig <- all[(all$padj < padj & (all$log2FoldChange>logFoldChange | all$log2FoldChange < (-logFoldChange))),]
diffSigout=diffSig[,colnames(diffSig)]
diffSigout=rbind(ID=colnames(diffSig),diffSig)
write.table(diffSigout,file="2.diff_mRNA.txt",sep="\t",quote=F,col.names=F)
#Write expression level of diff gene.
hmExp=outData[rownames(diffSig),]
diffGeneExpOut=rbind(id=paste0(colnames(hmExp),"_",Type),hmExp)
write.table(diffGeneExpOut, file="3.heatmap_mRNA.txt", sep="\t", quote=F, col.names=F)


#Output differential gene expression values.
diffGeneExp=data[row.names(diffSig),]
diffGeneExpOut=rbind(id=paste0(colnames(diffGeneExp),"_",Type), diffGeneExp)
write.table(diffGeneExpOut, file="diffGeneExp.txt", sep="\t", quote=F, col.names=F)

#Convert the adjusted p-values (adj.P.Val column) of differential genes to log10 scale.
all$logP <- -log10(all$padj)
#Add a new column named "Group".
all$Group = "not-significant"
#Set genes with adjusted p-values (adj.P.Val) less than 0.05 and log fold change (logFC) greater than 2 as significantly upregulated genes.
#Set genes with adjusted p-values (adj.P.Val) less than 0.05 and log fold change (logFC) less than 2 as significantly downregulated genes.
all$Group[which((all$padj <0.05)&(all$log2FoldChange > logFoldChange))]= "up-regulated"
all$Group[which((all$padj <0.05)&(all$log2FoldChange <(-logFoldChange)))]= "down-regulated"
#Check the number of upregulated and downregulated genes.
table(all$Group)
#Generate a new volcano plot.
#Add a new column named "Label".
all$Label = ""
#Sort the p-values of differentially expressed genes in ascending order.
all$Symbol=row.names(all)
all <- all[order(all$padj),]
#Among the highly expressed genes, select the 10 genes with the lowest adj.P.Val.
up.genes <- head(all$Symbol[which(all$Group == "up-regulated")],10)
#Among the lowly expressed genes, select the 10 genes with the lowest adj.P.Val.
down.genes <- head(all$Symbol[which(all$Group == "down-regulated")],10)
#Merge the upregulated genes and downregulated genes, and add them to the "Label" column.
deg.top10.genes <- c(as.character(up.genes),as.character(down.genes))
all$Label[match(deg.top10.genes,all$Symbol)] <- deg.top10.genes
#Change the colors of the volcano plot data points and adjust the axis labels to make the image more visually appealing.
pdf(file="mrnavol1.pdf",height=8,width=10)
xMax=max(abs(all$log2FoldChange))
ggscatter(all,x = "log2FoldChange",y = "logP",
          color = "Group",
          palette = c( "#2f5688","#BBBBBB","#CC0000"),size = 1,
          label = all$Label,font.label = 14,
          repel = T,
          xlim=c(-xMax,xMax),
          xlab = "log2FoldChange" ,
          ylab = "-log10(adj.P.value)",) + theme_base() +
  geom_hline(yintercept = 1.30,linetype="dashed") +
  geom_vline(xintercept = c(-logFoldChange,logFoldChange),linetype="dashed")
dev.off()

#Plotting a heatmap of differentially expressed genes.
#Write expression level of diff gene.
library(pheatmap)
hmExp=data[rownames(diffSig),]
#hmExp=expr[rownames(diffSig),]
Type=c(rep("Normal",conNum),rep("Tumor",treatNum))
names(Type)=colnames(diffExp)
Type=as.data.frame(Type)
pdf(file="mrnaHeatmap.pdf",height=12,width=15)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         show_rownames = F,
         scale="row",
         fontsize = 10,
         fontsize_row=4,
         fontsize_col=10)
dev.off()

