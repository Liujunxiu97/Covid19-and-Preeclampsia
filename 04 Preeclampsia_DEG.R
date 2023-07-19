#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#install.packages("pheatmap")
library(limma)      #Importing packages.
library(pheatmap)

inputFile="PreeclampsiaMatrix.txt"       #Input file.
logFCfilter=0.3           #logFC filtering threshold.
P.Valuefilter=0.05   #Adjusted p-value threshold after correction.
setwd("xxx")      #Set the working directory.

#Read the input file and organize the data from the input file.
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
rt=data

#Take the logarithm base 2 of the large values.
qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  rt[rt<0]=0
  rt=log2(rt+1)}
data=normalizeBetweenArrays(rt)

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

#Perform differential analysis.
Type=c(rep("con",conNum),rep("treat",treatNum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con","treat")
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
allDiff=topTable(fit2,adjust='fdr',number=200000)
allDiffOut=rbind(id=colnames(allDiff),allDiff)
write.table(allDiffOut, file="all.txt", sep="\t", quote=F, col.names=F)

#Output the corrected expression values.
outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData, file="normalize.txt", sep="\t", quote=F, col.names=F)

#Output the differential results.
diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & P.Value < P.Valuefilter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut, file="diff.txt", sep="\t", quote=F, col.names=F)

#Output the differential gene expression values.
diffGeneExp=data[row.names(diffSig),]
diffGeneExpOut=rbind(id=paste0(colnames(diffGeneExp),"_",Type), diffGeneExp)
write.table(diffGeneExpOut, file="diffGeneExp.txt", sep="\t", quote=F, col.names=F)

#Plotting a heatmap of differentially expressed genes.
geneNum=1000
diffSig=diffSig[order(as.numeric(as.vector(diffSig$logFC))),]
diffGeneName=as.vector(rownames(diffSig))
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(2*geneNum)){
    hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
    hmGene=diffGeneName
}
hmExp=data[hmGene,]
Type=c(rep("Con",conNum),rep("Treat",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf", width=10, height=8)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=7,
         fontsize_col=8)
dev.off()

#Performing a log10 transformation on the adjusted p-values (adj.P.Val column) of differentially expressed genes.
allDiff$logP <- -log10(allDiff$P.Value)
#Add a new column named "Group".
allDiff$Group = "not-significant"
#Set genes with adjusted p-values (adj.P.Val) less than 0.05 and log fold change (logFC) greater than 2 as significantly upregulated genes.
#Set genes with adjusted p-values (adj.P.Val) less than 0.05 and log fold change (logFC) less than 2 as significantly downregulated genes.
allDiff$Group[which((allDiff$P.Value <0.05)&(allDiff$logFC > logFCfilter))]= "up-regulated"
allDiff$Group[which((allDiff$P.Value <0.05)&(allDiff$logFC <(-logFCfilter)))]= "down-regulated"
#Check the number of upregulated and downregulated genes.
table(allDiff$Group)
#Generate a new volcano plot.
#Add a new column named "Label".
allDiff$Label = ""
#Sort the p-values of differentially expressed genes in ascending order.
allDiff$Symbol=row.names(allDiff)
allDiff <- allDiff[order(allDiff$P.Value),]
#Among the highly expressed genes, select the 10 genes with the lowest adj.P.Val.
up.genes <- head(allDiff$Symbol[which(allDiff$Group == "up-regulated")],10)
#Among the lowly expressed genes, select the 10 genes with the lowest adj.P.Val.
down.genes <- head(allDiff$Symbol[which(allDiff$Group == "down-regulated")],10)
#Merge the upregulated genes and downregulated genes and add them to the "Label" column.
deg.top10.genes <- c(as.character(up.genes),as.character(down.genes))
allDiff$Label[match(deg.top10.genes,allDiff$Symbol)] <- deg.top10.genes
#Change the colors of the volcano plot data points and adjust the axis labels to make the image more visually appealing.
pdf(file="mrnavol1.pdf",height=8,width=10)
xMax=max(abs(allDiff$logFC))
ggscatter(allDiff,x = "logFC",y = "logP",
          color = "Group",
          palette = c( "#2f5688","#BBBBBB","#CC0000"),size = 1,
          label = allDiff$Label,font.label = 14,
          repel = T,
          xlim=c(-xMax,xMax),
          xlab = "log2FoldChange" ,
          ylab = "-log10(adj.P.value)",) + theme_base() +
  geom_hline(yintercept = 1.30,linetype="dashed") +
  geom_vline(xintercept = c(-logFCfilter,logFCfilter),linetype="dashed")
dev.off()
dev.off()


