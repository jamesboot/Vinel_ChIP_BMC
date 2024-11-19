# Pathway analysis

# working directorly 
setwd("/Users/weiweijin/OneDrive - Queen Mary, University of London/DiffBind/pathway")
# setwd(system.file('extra',package='DiffBind'))

# clear workspace
rm(list = ls())

# Load packages
library(ChIPpeakAnno)
library(biomaRt) #to retrieve the gene names and properties
library(org.Hs.eg.db)
library(reactome.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library("tibble")
library(GenomicRanges)

# in/out parameters
dir_in <- 'input'
dir_out<- 'output_switched'

file_name <- 'ALL_K27me3-K27ac'

# readin files
gic_gain <- as.data.frame(read.table(file.path(dir_in, "ALL_K27ac_me3_switch_K27me3_ALL_v2_gic_loss.bed_K27ac_ALL_gic_gain.bed"),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
gic_gain <- gic_gain[,2:5]
gic_gain <- add_column(gic_gain, gic_gain[,3]-gic_gain[,2], .after = 3) 
colnames(gic_gain) <- c("chr","start","end","length","strand")

# pathway analysis etc
# data datapreprocessing
ensembl<-useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "http://www.ensembl.org",
                 ensemblRedirect = FALSE) # in case local host is down
annaData<-getAnnotation(ensembl, featureType = "ExonPlusUtr")
granges_gain <- GRanges(gic_gain)
gain.peaksAnno=annotatePeakInBatch(granges_gain, AnnotationData=annaData, output="overlapping", bindingType = "fullRange", bindingRegion = c(-5000, 1))

# TE
transposableElements<-read.table("repeats.bed")
colnames(transposableElements)<-c("chr","start","end","feature","width", "strand")
transposableElements$chr<-gsub("chr","",as.character(transposableElements$chr))
transposableElements<-makeGRangesFromDataFrame(transposableElements,keep.extra.columns=TRUE)
gain.hits <- findOverlaps(transposableElements,gain.peaksAnno)
gain.idx <- unique(subjectHits(gain.hits))
gain.TE <- gain.peaksAnno[gain.idx]

values_gain<- gain.TE$ensembl_gene_id
list.attr.ens=c("ensembl_gene_id","external_gene_name")
geneID_gain<-getBM(attributes=list.attr.ens, filters="ensembl_gene_id", values=values_gain, mart=mart)
matching_names_gain=vector()
for (ln in seq(length(gain.TE))) {
  matching_names_gain[ln]<-geneID_gain$external_gene_name[match(values_gain[ln],geneID_gain$ensembl_gene_id)]
}
data_gain.geneName <- gain.TE
data_gain.geneName$geneName <- matching_names_gain
write.table(data_gain.geneName, file=file.path(dir_out,paste(file_name,"_gain_withGeneName_TE.csv", sep = "")), sep = ",", quote = FALSE, row.names=F)

# get gene names and save
values_gain<- gain.peaksAnno$ensembl_gene_id
list.attr.ens=c("ensembl_gene_id","external_gene_name")
geneID_gain<-getBM(attributes=list.attr.ens, filters="ensembl_gene_id", values=values_gain, mart=ensembl)
matching_names_gain=vector()
for (ln in seq(length(gain.peaksAnno))) {
  matching_names_gain[ln]<-geneID_gain$external_gene_name[match(values_gain[ln],geneID_gain$ensembl_gene_id)]
}
data_gain.geneName <- gain.peaksAnno
data_gain.geneName$geneName <- matching_names_gain
write.table(data_gain.geneName, file=file.path(dir_out,paste(file_name,"_gain_withGeneName.csv", sep = "")), sep = ",", quote = FALSE, row.names=F)

# Peak distribution over genomic features
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
gain.featuresDist<-assignChromosomeRegion(granges_gain, nucleotideLevel=FALSE, precedence=c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs","Exons", "Introns"), TxDb=txdb)

pdf(file=file.path(dir_out,paste(file_name,"_peaks_featuresDistr_gic_gain.pdf", sep = "")), height=8, width=8, compress=TRUE)
par(mar=c(5, 10, 4, 2) + 0.1)
barplot(gain.featuresDist$percentage, las=1, horiz=T)
dev.off()

# GO ontologies
gain.go <- getEnrichedGO(gain.peaksAnno, orgAnn="org.Hs.eg.db", maxP=.1, minGOterm=10, multiAdjMethod="BH", condense=TRUE)

# Save GO ontologies results
write.table(gain.go$bp, file=file.path(dir_out,paste(file_name,"_go_bp_gain_gic.csv", sep = "")), sep=",", quote=F, row.names=F, col.names=T)
write.table(gain.go$mf, file=file.path(dir_out,paste(file_name,"_go_mf_gain_gic.csv", sep = "")), sep=",", quote=F, row.names=F, col.names=T)
write.table(gain.go$cc, file=file.path(dir_out,paste(file_name,"_go_cc_gain_gic.csv", sep = "")), sep=",", quote=F, row.names=F, col.names=T)

# REACTOME pathways
gain.pathways <- getEnrichedPATH(gain.peaksAnno, "org.Hs.eg.db", "reactome.db", maxP=.05)

# REACTOME pathways: save data
write.table(gain.pathways, file=file.path(dir_out,paste(file_name,"_reactome_pathyway_gain_gic.csv", sep = "")), sep=",", quote=F, row.names=F, col.names=T)

