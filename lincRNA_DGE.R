
### First, clear the workspace to free up memory

rm(list=ls())
options(expressions=500000)
### if the packages have not been previously loaded remove the #'s before the following lines

source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("simpleaffy")
biocLite("limma")
biocLite("hgu133plus2.db")
biocLite("annotate")
biocLite("genefilter")

### Library() tells R to load the following packages for our purposes.

library(genefilter)
library(hgu133plus2.db)
library(annotate)
library(simpleaffy)
library(limma)
library(gplots)
library(corrplot)
library(dendextend)

### Read in a 3 column text file conaining the path to the files, the sample name used as the column header, and the ALL target group to which they belong
### Then store the factors for later use

pheno<-read.table("/nobackup/shilab/Data/StJude/Leukemia_Subtypes/phenotypes/pheno_BCR_Phlike_crlf2.txt",sep="\t",header=TRUE)
samples<-as.factor(pheno$Target)

### The simpleaffy package will use the paths from the text file to find the microarray CEL files, normalize them using Robust multi-array average.
### We will use the factors created earlier to reassign the empty "mol.biol" list for the newly created eset data structure,
### this will allow us to produce a color bar along the heatmap reflecting which group is which.

celfiles<-read.affy(covdesc="/phenotypes/pheno_BCR_Phlike_crlf2.txt")
celfiles.rma<-rma(celfiles)
celfiles.rma$mol.biol<-samples
eset<-exprs(celfiles.rma)
head(rownames(eset))

celfiles.filtered<-nsFilter(celfiles.rma, require.extrez=FALSE, remove.dupEntrez= FALSE)
design<-model.matrix(~0+samples)
colnames(design)<-levels(samples)
fit<-lmFit(exprs(celfiles.filtered$eset), design)

### The contrast phenotypes should be written exactly as the in the Targets.
contrast.matrix<-makeContrasts(diff="Ph_like_CRLF2-BCR_ABL1",levels=design)
fit2<-contrasts.fit(fit,contrast.matrix)
eBfit<-eBayes(fit2)

### Here we create a function to differentiate the samples using hex-codes based on the mol.biol data from the eset.

color.map<-function(mol.biol){ if (mol.biol=="Ph_like_CRLF2") {"#FFE135" } else if (mol.biol=="Ph_like_non_CRLF2"){ "#007FFF" } else {"#FF00FF"} }

### The hgu133plus2 affy chip contains over 54,000 probes. 
### topTable provides a table of just the top n probes, in this case 40,000, which are fdr corrected for multiple testing error and annotated.

mytopTable1<-topTable(eBfit,number=Inf, adjust="fdr")
gene.symbols<-getSYMBOL(row.names(mytopTable1),"hgu133plus2")
results<-cbind(gene.symbols,mytopTable1)

### The probes are now collapsed to give a gene-level analysis of differential expression
### The head command allows us to preview what each step looks like as we go and is written in the output file *.Rout
### This is good practice for debugging.
### Finally we use a two-fold Beta value and 1% of a type-one error as our selection criteria

myresmax<-findLargest(rownames(results),results$AveExpr,"hgu133plus2")
head(myresmax)
selection<-results[myresmax,]
head(selection)
selected<-subset(selection,selection$B >= 1.5 & selection$adj.P.Val < 0.01)

### Let's see how many genes fit our criteria and will be included in the heatmap

dim(selected)
head(selected)
selected<-selected[with(selected ,order(-abs(selected$B))), ]
head(selected)
selected<-selected[1:50,]
dim(selected)
head(selected)
esetSel<-celfiles.filtered$eset[rownames(selected), ]
dim(esetSel)


write.table(eset,file="/users/cnodzak/esetBCR_PhlikeCRLF2.txt", row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t")

#write.table(selection,file="/nobackup/shilab/Data/StJude/Leukemia_Subtypes/esets/selection_BCR.txt", row.names=TRUE, col.names=TRUE,quote=FALSE,sep="\t")
