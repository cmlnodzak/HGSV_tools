


crlf2<-read.table("crlf2_genes.bed.intersect.txt.vcf.hetSNPs.vcf",sep="\t")
counts<-read.table("crlf2_genes.bed.intersect.txt.vcf.hetSNPs.vcf.counts.txt",sep="\t")

mymerge<-merge(crlf2,counts, by.x="pos",by.y="POS")


mymerge$pval<-apply(mymerge,1,function(x) binom.test(as.numeric(x[4]),as.numeric(x[7]))$p.value)

pval<-list()
for (i in 1:nrow(mymerge)){
  if (mymerge[i,5]=="A") { 
    if(mymerge[i,7]==0)
      {pval[i]<-"NA"}
    else{
    pval[i]<-binom.test(as.numeric(mymerge[i,7]),as.numeric(mymerge[i,6]))$p.value
  }
  } else if (mymerge[i,5]=="C"){
    if(mymerge[i,8]==0)
      {pval[i]<-"NA"}
    else{
      pval[i]<-binom.test(as.numeric(mymerge[i,8]),as.numeric(mymerge[i,6]))$p.value
    }
  } else if (mymerge[i,5]=="G"){
    if(mymerge[i,9]==0)
      {pval[i]<-"NA"}
    else{
      pval[i]<-binom.test(as.numeric(mymerge[i,9]),as.numeric(mymerge[i,6]))$p.value
    }
  } else if (mymerge[i,5]=="T"){
    if(mymerge[i,10]==0)
      {pval[i]<-"NA"}
    else{
      pval[i]<-binom.test(as.numeric(mymerge[i,10]),as.numeric(mymerge[i,6]))$p.value
    }
  } else{
    pval[i]<-"NA"
  }
  }
pval<-unlist(pval)
mymerge<-cbind.data.frame(mymerge,pval)
mymerge$FDR<-p.adjust(pval,method="fdr")
ASEsignif<-mymerge[mymerge$FDR<=0.01,]


write.table(ASEsignif,"crlf2.ase.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)


