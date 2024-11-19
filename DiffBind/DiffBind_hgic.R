# DiffBind

# working directorly 
setwd("/Users/weiweijin/OneDrive - Queen Mary, University of London/DiffBind")
# setwd(system.file('extra',package='DiffBind'))

# clear workspace
rm(list = ls())

# Load packages
library(DiffBind)
library(ChIPpeakAnno)
library(biomaRt) #to retrieve the gene names and properties
library(tidyverse)
library(dplyr)
library(org.Hs.eg.db)
library(reactome.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(parallel)
BPPARAM=BiocParallel::SerialParam()

# import data & setting
subgroup <- "RTKII"
antibody <- "K4me3"
file_name <- paste(antibody, "_", subgroup, sep = "")
input_dir <- "input"
output_dir <- "output"
data_input <- dba(sampleSheet=file.path(input_dir,paste(file_name,".csv", sep = "")))
pdf(file=file.path(output_dir, paste(file_name,"_org.pdf")), height=8, width=8, compress=TRUE)
plot(data_input)
dev.off();

# count reads
data_count <- dba.count(data_input, summits = 1000)
info <- dba.show(data_count)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP, PeakReads=round(info$Reads * info$FRiP))
rownames(libsizes) <- info$ID
pdf(file=file.path(output_dir,paste(file_name,"_count.pdf", sep = "")), height=8, width=8, compress=TRUE)
plot(data_count)
dev.off();

# Normalizing the data
data_norm <- dba.normalize(data_count, bRetrieve=TRUE)

# Establishing a model design and contrast
data_design <- dba.contrast(data_count, minMembers=2, reorderMeta=list(Condition="NSC"))

# performing differntial analysis
data_analysed <- dba.analyze(data_design, method=DBA_ALL_METHODS)
dba.show(data_analysed, bContrasts=TRUE)
pdf(file=file.path(output_dir,paste(file_name,"_analysed.pdf", sep = "")), height=8, width=8, compress=TRUE)
plot(data_analysed, contrast=1)
dev.off();

# retrieve the differentially bound sites
data.DB <- dba.report(data_analysed, contrast=1)
sum(data.DB$Fold>0)
sum(data.DB$Fold<0)
write.table(data.DB, file = file.path(output_dir,paste(file_name,"_sum_table.csv", sep = "")), sep = ",", quote = FALSE, row.names = F)

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

# save data to bed files
out <- as.data.frame(data.DB)
gic_gain <- out %>% 
  dplyr::filter(FDR < 0.05 & Fold > 0) %>% 
  dplyr::select(seqnames, start, end)
write.table(gic_gain, file=file.path(output_dir,paste(file_name,"_gic_gain.bed", sep = "")), sep="\t", quote=F, row.names=F, col.names=T)
gic_loss <- out %>% 
  dplyr::filter(FDR < 0.05 & Fold < 0) %>% 
  dplyr::select(seqnames, start, end)
write.table(gic_loss, file=file.path(output_dir,paste(file_name,"_gic_loss.bed", sep = "")), sep="\t", quote=F, row.names=F, col.names=T)

# Enriched GO / REACTOME terms
# Peak distribution over genomic features
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peaks.featuresDist<-assignChromosomeRegion(data.peaksAnno, nucleotideLevel=FALSE, precedence=c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs","Exons", "Introns"), TxDb=txdb)

pdf(file=file.path(output_dir,paste(file_name,"_peaks_featuresDistr.pdf", sep = "")), height=8, width=8, compress=TRUE)
par(mar=c(5, 10, 4, 2) + 0.1)
barplot(peaks.featuresDist$percentage, las=1, horiz=T)
dev.off()

# GO ontologies
# peaks.go <- getEnrichedGO(data.peaksAnno, orgAnn="org.Hs.eg.db", maxP=.1, minGOterm=10, multiAdjMethod="BH", condense=TRUE)
# 
# # Save GO ontologies results
# write.table(peaks.go$bp, file=file.path(output_dir,paste(file_name,"_go_bp.csv", sep = "")), sep=",", quote=F, row.names=F, col.names=T)
# write.table(peaks.go$mf, file=file.path(output_dir,paste(file_name,"_go_mf.csv", sep = "")), sep=",", quote=F, row.names=F, col.names=T)
# write.table(peaks.go$cc, file=file.path(output_dir,paste(file_name,"_go_cc.csv", sep = "")), sep=",", quote=F, row.names=F, col.names=T)

# REACTOME pathways
peaks.pathways <- getEnrichedPATH(data.peaksAnno, "org.Hs.eg.db", "reactome.db", maxP=.05)

# REACTOME pathways: save data
write.table(peaks.pathways, file=file.path(output_dir,paste(file_name,"_reactome_pathyway.csv", sep = "")), sep=",", quote=F, row.names=F, col.names=T)

# Venn diagrams
pdf(file=file.path(output_dir,paste(file_name,"_venn_gain_loss.pdf", sep = "")), height=8, width=8, compress=TRUE)
dba.plotVenn(data_analysed, contrast=1, bDB=TRUE, bGain=TRUE, bLoss=TRUE, bAll=FALSE)
dev.off();

pdf(file=file.path(output_dir,paste(file_name,"_venn_methods.pdf", sep = "")), height=8, width=8, compress=TRUE)
dba.plotVenn(data_analysed, contrast=1, method=DBA_ALL_METHODS)
dev.off();

# PCA
pdf(file=file.path(output_dir,paste(file_name,"_pca_pre.pdf", sep = "")), height=8, width=8, compress=TRUE)
dba.plotPCA(data_analysed,DBA_TISSUE,label=DBA_CONDITION)
dev.off();

pdf(file=file.path(output_dir,paste(file_name,"_pca_post.pdf", sep = "")), height=8, width=8, compress=TRUE)
dba.plotPCA(data_analysed, contrast=1, label=DBA_TISSUE)
dev.off();

# MA plots
pdf(file=file.path(output_dir,paste(file_name,"_ma.pdf", sep = "")), height=8, width=8, compress=TRUE)
dba.plotMA(data_analysed)
dev.off();

# Vocano plots
pdf(file=file.path(output_dir,paste(file_name,"_vocano.pdf", sep = "")), height=8, width=8, compress=TRUE)
dba.plotVolcano(data_analysed)
dev.off();

# Boxplots
pdf(file=file.path(output_dir,paste(file_name,"_boxplot.pdf", sep = "")), height=8, width=8, compress=TRUE)
pvals <- dba.plotBox(data_analysed)
dev.off();

# Heatmaps
pdf(file=file.path(output_dir,paste(file_name,"_corvals.pdf", sep = "")), height=8, width=8, compress=TRUE)
corvals <- dba.plotHeatmap(data_analysed)
dev.off();

hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)
pdf(file=file.path(output_dir,paste(file_name,"_readscore.pdf", sep = "")), height=8, width=8, compress=TRUE)
readscores <- dba.plotHeatmap(data_analysed, contrast=1, correlations=FALSE, scale="row", colScheme = hmap)
dev.off();

# Profiling and Profile Heatmaps
all_profiles <- dba.plotProfile(data_analysed)
pdf(file=file.path(output_dir,paste(file_name,"_profile_sample.pdf", sep = "")), height=8, width=8, compress=TRUE)
dba.plotProfile(all_profiles)
dev.off();

merge_profiles <- dba.plotProfile(data_analysed, merge = c(DBA_TISSUE,DBA_REPLICATE))
pdf(file=file.path(output_dir,paste(file_name,"_profile_merge.pdf", sep = "")), height=8, width=8, compress=TRUE)
dba.plotProfile(merge_profiles)
dev.off();

# save data
save.image(file = file.path(output_dir,paste(file_name,".RData", sep = "")))

