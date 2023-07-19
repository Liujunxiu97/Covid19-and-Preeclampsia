#install.packages("pROC")
library(pROC)                  #Import packages.
expFile="3.heatmap_mRNA.txt"      #Expression data file.
geneFile="hub 10.txt"      #Intersection gene list file.
setwd("xxx")    #Set the working directory.
#Read the input file and organize the data from the input file.
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
y=gsub("(.*)\\_(.*)", "\\2", colnames(rt))
y=ifelse(y=="con", 0, 1)

#Read the gene list file.
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)

#Loop through the intersection genes and plot ROC curves.
for(x in as.vector(geneRT[,1])){
#Plot the ROC curve.
	roc1=roc(y, as.numeric(rt[x,]))
	ci1=ci.auc(roc1, method="bootstrap")
	ciVec=as.numeric(ci1)
	pdf(file=paste0("ROC.",x,".pdf"), width=5, height=5)
	plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=x)
	text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
	dev.off()
}

