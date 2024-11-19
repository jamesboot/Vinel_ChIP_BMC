# peak annotation 

# working directorly 
setwd("/Users/weiweijin/OneDrive - Queen Mary, University of London/DiffBind/peak_annotation")
# setwd(system.file('extra',package='DiffBind'))

# clear workspace
rm(list = ls())

# Load packages
library(ChIPpeakAnno)
library(biomaRt) #to retrieve the gene names and properties


# import data & setting



# annotate the differentially bound sites
data.peaksAnno=annotatePeakInBatch(data.DB, AnnotationData=TSS.human.GRCh38)
write.table(data.peaksAnno, file=file.path(output_dir,paste(file_name,"_anno.csv", sep = "")), sep = ",", quote = FALSE, row.names=F)

# get gene names
ensembl<-useMart("ensembl", dataset="hsapiens_gene_ensembl") 
values<- data.peaksAnno$feature
list.attr.ens=c("ensembl_gene_id","external_gene_name")
geneID<-getBM(attributes=list.attr.ens, filters="ensembl_gene_id", values=values, mart=ensembl)
matching_names=vector()
for (ln in seq(length(data.peaksAnno))) {
  matching_names[ln]<-geneID$external_gene_name[match(values[ln],geneID$ensembl_gene_id)]
}
data.geneName <- data.peaksAnno
data.geneName$geneName <- matching_names
write.table(data.geneName, file=file.path(output_dir,paste(file_name,"_withGeneName.csv", sep = "")), sep = ",", quote = FALSE, row.names=F)
