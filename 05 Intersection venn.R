#install.packages("VennDiagram")
library(VennDiagram)        #Importing packages.
setwd("xxx")    #Set the working directory.
files=dir()     #Get all files in the directory.
files=grep("txt",files,value=T)     #Extract files with the .txt extension.
geneList=list()

#Read the gene information from all .txt files and save it to geneList.
for(i in 1:length(files)){
    inputFile=files[i]
	if(inputFile=="intersect.txt"){next}
    rt=read.table(inputFile,header=F)
    header=unlist(strsplit(inputFile,"\\.|\\-"))
    geneList[[header[1]]]=as.vector(rt[,1])
    uniqLength=length(unique(as.vector(rt[,1])))
    print(paste(header[1],uniqLength,sep=" "))
}

#Plot a Venn diagram.
venn.plot=venn.diagram(geneList,filename=NULL,fill=rainbow(length(geneList)))
pdf(file="venn.pdf",width=15,height=15)
grid.draw(venn.plot)
dev.off()

#Save the intersecting genes.
intersectGenes=Reduce(intersect,geneList)
write.table(file="intersect.txt",intersectGenes,sep="\t",quote=F,col.names=F,row.names=F)

