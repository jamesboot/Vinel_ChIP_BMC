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
dir_out<- 'output_potential_enhancer'

file_name <- 'MES_K27ac-K4me3'

# readin files
gic_gain <- as.data.frame(read.table(file.path(dir_in, "MES_K27ac_K4_K27ac_MES_gic_gain.bed"),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
gic_gain <- gic_gain[,2:5]
gic_gain <- add_column(gic_gain, gic_gain[,3]-gic_gain[,2], .after = 3) 
colnames(gic_gain) <- c("chr","start","end","length","strand")

loss_gain <- as.data.frame(read.table(file.path(dir_in, "MES_K27ac_K4_K27ac_MES_gic_loss.bed"),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
loss_gain <- loss_gain[,2:5]
loss_gain <- add_column(loss_gain, loss_gain[,3]-loss_gain[,2], .after = 3)
colnames(loss_gain) <- c("chr","start","end","length","strand")

# pathway analysis etc
# data datapreprocessing
ensembl<-useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "http://www.ensembl.org",
                 ensemblRedirect = FALSE) # in case local host is down
annaData_all <-getAnnotation(ensembl, featureType = "ExonPlusUtr")
annaData_transcript <-getAnnotation(ensembl, featureType = "transcript")
annaData_exon <-getAnnotation(ensembl, featureType = "Exon")
granges_gain <- GRanges(gic_gain)

gain.peaksAnno_all= annotatePeakInBatch(granges_gain, AnnotationData=annaData_all, 
                                        output="overlapping", bindingType = "fullRange", bindingRegion = c(-5000, 1))
gain.peaksAnno_all$feature <- NULL
gain.peaksAnno_all$feature.ranges <- NULL
gain.peaksAnno_all$distance <- NULL
gain.peaksAnno_all$insideFeature <- NULL
gain.peaksAnno_all$ensembl_exon_id <- NULL
gain.peaksAnno_all$utr5start <- NULL
gain.peaksAnno_all$utr5end <- NULL
gain.peaksAnno_all$utr3start <- NULL
gain.peaksAnno_all$utr3end <- NULL
gain.peaksAnno_all$distanceToSite <- NULL
gain.peaksAnno_all$description <- NULL
hits<-findOverlaps(granges_gain,gain.peaksAnno_all)
idx <- unique(queryHits(hits))
gain.peaksAnno_intergenic <- annotatePeakInBatch(granges_gain[-idx], AnnotationData=TSS.human.GRCh38)

gain.peaksAnno_transcript = annotatePeakInBatch(granges_gain, AnnotationData = annaData_transcript, 
                                                output="overlapping", FeatureLocForDistance="geneEnd", bindingRegion = c(-1, 1))
gain.peaksAnno_exon = annotatePeakInBatch(granges_gain, AnnotationData = annaData_exon, 
                                                output="overlapping", FeatureLocForDistance="geneEnd", bindingRegion = c(-1, 1))
gain.peaksAnno_exon$feature <- NULL
gain.peaksAnno_exon$feature.ranges <- NULL
gain.peaksAnno_exon$insideFeature <- NULL
gain.peaksAnno_exon$description <- NULL
hits2<-findOverlaps(gain.peaksAnno_transcript,gain.peaksAnno_exon)
idx2 <- unique(queryHits(hits2))
gain.peaksAnno_intron <- gain.peaksAnno_transcript[-idx2]

granges_loss <- GRanges(loss_gain)
loss.peaksAnno_all=annotatePeakInBatch(granges_loss, AnnotationData=annaData_all,
                                       output="overlapping", bindingType = "fullRange", bindingRegion = c(-5000, 1))
hits<-findOverlaps(granges_loss,loss.peaksAnno_all)
idx <- unique(queryHits(hits))
loss.peaksAnno_intergenic <- annotatePeakInBatch(granges_loss[-idx], AnnotationData=TSS.human.GRCh38)

loss.peaksAnno_transcript = annotatePeakInBatch(granges_loss, AnnotationData = annaData_transcript, 
                                                output="overlapping", FeatureLocForDistance="geneEnd", bindingRegion = c(-1, 1))
loss.peaksAnno_exon = annotatePeakInBatch(granges_loss, AnnotationData = annaData_exon, 
                                          output="overlapping", FeatureLocForDistance="geneEnd", bindingRegion = c(-1, 1))
hits2<-findOverlaps(loss.peaksAnno_transcript,loss.peaksAnno_exon)
idx2 <- unique(queryHits(hits2))
loss.peaksAnno_intron <- loss.peaksAnno_transcript[-idx2]

# get gene names and save
values_gain<- gain.peaksAnno_intergenic$feature
values_loss<- loss.peaksAnno_intergenic$feature
list.attr.ens=c("ensembl_gene_id","external_gene_name")
geneID_gain<-getBM(attributes=list.attr.ens, filters="ensembl_gene_id", values=values_gain, mart=ensembl)
geneID_loss<-getBM(attributes=list.attr.ens, filters="ensembl_gene_id", values=values_loss, mart=ensembl)
matching_names_gain=vector()
matching_names_loss=vector()
for (ln in seq(length(gain.peaksAnno_intergenic))) {
  matching_names_gain[ln]<-geneID_gain$external_gene_name[match(values_gain[ln],geneID_gain$ensembl_gene_id)]
}
for (ln in seq(length(loss.peaksAnno_intergenic))) {
  matching_names_loss[ln]<-geneID_loss$external_gene_name[match(values_loss[ln],geneID_loss$ensembl_gene_id)]
}
data_gain_intergenic.geneName <- gain.peaksAnno_intergenic
data_loss_intergenic.geneName <- loss.peaksAnno_intergenic
data_gain_intergenic.geneName$geneName <- matching_names_gain
data_loss_intergenic.geneName$geneName <- matching_names_loss
write.table(data_gain_intergenic.geneName, file=file.path(dir_out,paste(file_name,"_gain_withGeneName.csv", sep = "")), sep = ",", quote = FALSE, row.names=F)
write.table(data_loss_intergenic.geneName, file=file.path(dir_out,paste(file_name,"_loss_withGeneName.csv", sep = "")), sep = ",", quote = FALSE, row.names=F)

values_gain<- gain.peaksAnno_intron$ensembl_gene_id
values_loss<- loss.peaksAnno_intron$ensembl_gene_id
geneID_gain<-getBM(attributes=list.attr.ens, filters="ensembl_gene_id", values=values_gain, mart=ensembl)
geneID_loss<-getBM(attributes=list.attr.ens, filters="ensembl_gene_id", values=values_loss, mart=ensembl)
matching_names_gain=vector()
matching_names_loss=vector()
for (ln in seq(length(gain.peaksAnno_intron))) {
  matching_names_gain[ln]<-geneID_gain$external_gene_name[match(values_gain[ln],geneID_gain$ensembl_gene_id)]
}
for (ln in seq(length(loss.peaksAnno_intron))) {
  matching_names_loss[ln]<-geneID_loss$external_gene_name[match(values_loss[ln],geneID_loss$ensembl_gene_id)]
}
data_gain_intron.geneName <- gain.peaksAnno_intron
data_loss_intron.geneName <- loss.peaksAnno_intron
data_gain_intron.geneName$geneName <- matching_names_gain
data_loss_intron.geneName$geneName <- matching_names_loss
write.table(data_gain_intron.geneName, file=file.path(dir_out,paste(file_name,"_gain_withGeneName.csv", sep = "")), sep = ",", quote = FALSE, row.names=F)
write.table(data_loss_intron.geneName, file=file.path(dir_out,paste(file_name,"_loss_withGeneName.csv", sep = "")), sep = ",", quote = FALSE, row.names=F)


# Peak distribution over genomic features
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
gain.featuresDist<-assignChromosomeRegion(granges_gain, nucleotideLevel=FALSE, precedence=c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs","Exons", "Introns"), TxDb=txdb)

pdf(file=file.path(dir_out,paste(file_name,"_peaks_featuresDistr_gic_gain.pdf", sep = "")), height=8, width=8, compress=TRUE)
par(mar=c(5, 10, 4, 2) + 0.1)
barplot(gain.featuresDist$percentage, las=1, horiz=T)
dev.off()

loss.featuresDist<-assignChromosomeRegion(granges_loss, nucleotideLevel=FALSE, precedence=c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs","Exons", "Introns"), TxDb=txdb)

pdf(file=file.path(dir_out,paste(file_name,"_peaks_featuresDistr_gic_loss.pdf", sep = "")), height=8, width=8, compress=TRUE)
par(mar=c(5, 10, 4, 2) + 0.1)
barplot(loss.featuresDist$percentage, las=1, horiz=T)
dev.off()

# GO ontologies
gain.go <- getEnrichedGO(gain.peaksAnno, orgAnn="org.Hs.eg.db", maxP=.1, minGOterm=10, multiAdjMethod="BH", condense=TRUE)

loss.go <- getEnrichedGO(loss.peaksAnno, orgAnn="org.Hs.eg.db", maxP=.1, minGOterm=10, multiAdjMethod="BH", condense=TRUE)

# Save GO ontologies results
write.table(gain.go$bp, file=file.path(dir_out,paste(file_name,"_go_bp_gain_gic.csv", sep = "")), sep=",", quote=F, row.names=F, col.names=T)
write.table(gain.go$mf, file=file.path(dir_out,paste(file_name,"_go_mf_gain_gic.csv", sep = "")), sep=",", quote=F, row.names=F, col.names=T)
write.table(gain.go$cc, file=file.path(dir_out,paste(file_name,"_go_cc_gain_gic.csv", sep = "")), sep=",", quote=F, row.names=F, col.names=T)

write.table(loss.go$bp, file=file.path(dir_out,paste(file_name,"_go_bp_loss_gic.csv", sep = "")), sep=",", quote=F, row.names=F, col.names=T)
write.table(loss.go$mf, file=file.path(dir_out,paste(file_name,"_go_mf_loss_gic.csv", sep = "")), sep=",", quote=F, row.names=F, col.names=T)
write.table(gain.go$cc, file=file.path(dir_out,paste(file_name,"_go_cc_loss_gic.csv", sep = "")), sep=",", quote=F, row.names=F, col.names=T)

# REACTOME pathways
gain.pathways <- getEnrichedPATH(gain.peaksAnno, "org.Hs.eg.db", "reactome.db", maxP=.05)

loss.pathways <- getEnrichedPATH(loss.peaksAnno, "org.Hs.eg.db", "reactome.db", maxP=.05)

# REACTOME pathways: save data
write.table(gain.pathways, file=file.path(dir_out,paste(file_name,"_reactome_pathyway_gain_gic.csv", sep = "")), sep=",", quote=F, row.names=F, col.names=T)

write.table(loss.pathways, file=file.path(dir_out,paste(file_name,"_reactome_pathyway_loss_gic.csv", sep = "")), sep=",", quote=F, row.names=F, col.names=T)

