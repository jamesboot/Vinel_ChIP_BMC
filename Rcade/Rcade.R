# Rcade for intergation of ChIP seq and DE 

#to set the Working Directory; this depends on the local machine
setwd("/media/weiweijin/Stan/Rcade")

# clear workspace
rm(list = ls())

# Load packages
library(Rcade)
library(biomaRt)
library(parallel)
library(rgl)
library(tidyverse)

# find the directory of the sample data
# dir <- file.path(system.file("extdata", package="Rcade"), "STAT1")

# get the sample DE matrix 
# DE matrix: Bayesian information about the DE status of each gene. Any Bayesian source can be used for
# this purpose:
#   • limma: Full DE results obtained with topTable(..., number=Inf)
# • baySeq.
# At time of writing, edgeR and DEseq do not supply Bayesian output.
# The DE data must contain the following fields:
#   geneID - gene IDs used to link DE results to genes.
# logFC - The log fold change associated with each gene.
# B - B values (log-odds).
de_dir <- "de_input"
de_files <- list.files(path = de_dir,
                    pattern="*.csv",
                    full.names = TRUE)
print(de_files)

# Which col in DE matrix is important for Rcade
# This is done using
# an object of the form list(RcadeField1=DEdataField1, ...). You may omit the name RcadeField1
# if it is the same as DEdataField1. Any fields that are specified in addition to the three required fields
# above will not be manipulated in the analysis, but will be carried through and appear in the output.
DElookup <- list(GeneID="ENSG", logFC="logFC", B="B",
                 "Names")

# Annotation information
anno <- getBM(
  attributes= c("ensembl_gene_id", "chromosome_name",
                "transcript_start", "transcript_end", "strand"),
  mart= useDataset("hsapiens_gene_ensembl", useMart("ensembl")), useCache = F
)
##order, to reduce size of ChIPannoZones object later
anno <- anno[order(anno$chromosome_name),]
##use appropriate column names
colnames(anno) <- c("ENSG","chr","start","end","str")

# In the ChIP-seq analysis, Rcade performs its analysis based on the counts in user-defined bin regions. These
# regions are specified with a GRanges object from the GenomicRanges package.
# A common requirement is to define bins about genomic annotation features: Rcade provides simple functionality
# to generate an appropriate GRanges object, through the defineBins() function. Since STAT1 is a
# promoter-bound TF, we define bins about Ensembl-derived transcription start sites.
ChIPannoZones <- defineBins(anno, zone=c(-2000, 2000), geneID="ENSG")
# The zone=c(-1500,1500) argument defines the zone of interest: it starts 1500bp 50 of each TSS (-1500),
# and ends 1500bp 30 of each TSS (1500). This zone was chosen based on mapping peak-calls to TSSs using
# the ChIPpeakAnno package, and plotting the positional distribution of these peak-calls.

# Prior specification
# We specify DE.prior = 0.01 because that is the prior probability used by limma, the package that our DE
# analysis was performed with. The remaining settings ensure that genes with ChIP-seq signal have a higher
# than average probability of DE ("D|C" = 0.05), whereas genes without ChIP-seq signal have a lower that
# average probability of DE ("D|notC" = 0.005).
# (The values 0.05 and 0.005 were selected arbitrarily, before analysis. An advanced user might select prior
#   in a less arbitrary manner – for example, by looking at the overlap between ChIP-seq and DE in “similar"
# datasets. We do not go into further details of such an analysis here.)
# If you do not supply this information, you may get unreasonably small B values at the end of the analysis.
DE.prior = 0.01
prior.mode = "keepChIP"
prior = c("D|C" = 0.05, "D|notC" = 0.005)

# Analysis function 
# parallelse
# n_jobs <- as.integer(Sys.getenv("NSLOTS"))
n_jobs <- 12
print(paste("numer of cores: ", n_jobs))
cl <- makeCluster(n_jobs, "SOCK")

chip_dir <- "chip_input"

# i = as.integer(Sys.getenv("SGE_TASK_ID"))
for (i in 1:3){
print(paste("opening de file: ", de_files[i]))
DE <- read.csv(de_files[i])

output_dir <- sub("de_input/RNA_", "output/", de_files)

# ChIP data: .bam and .bai files: Sets of aligned reads. These reads should have already undergone sequence-level
# pre-processing, such as any read trimming and adaptor removal that may be required. Moreover,
# they should have appropriate index files - for example, using the indexBam function in the package
# Rbamtools.

# Targets information: A matrix containing information about the .Bam files to be used in the analysis.
# Required fields:
#   fileID – ID associated with the file.
# sampleID – ID associated with the sample that was sequenced. Technical replicates of the same
# population should have the same sampleID.
# factor – The antibody used in the experiment. Control files should be labelled "Input".
# filepath – The .bam file’s name/path.
# Optional fields:
#   shift – Half of the fragment length. That is, before counting, Rcade will shift reads on the positive
# strand forwards by shift, and reads on the negative strand backwards by shift.
target_files <- sub("de_input/RNA", "other_input/targets", de_files)
for (j in 1:4){
  if(j == 1){ChIPannoZones <- defineBins(anno, zone=c(-1500, 1500), geneID="ENSG")}
  else if(j ==2){ChIPannoZones <- defineBins(anno, zone=c(-4500, 4500), geneID="ENSG")}
  else if(j ==3){ChIPannoZones <- defineBins(anno, zone=c(-2000, 2000), geneID="ENSG")}
  else {ChIPannoZones <- defineBins(anno, zone=c(-1500, 1500), geneID="ENSG")}  
  target_file <- sub(".csv", paste("_", as.character(j), ".csv", sep = ""), target_files[i])
  print(paste("opening target file: ", target_file))
  targets <- read.csv(target_file, as.is = TRUE)

  Rcade <- RcadeAnalysis(DE, ChIPannoZones, annoZoneGeneidName="ENSG",
                        ChIPtargets=targets, ChIPfileDir = chip_dir,
                        cl=cl, DE.prior=DE.prior, prior.mode=prior.mode, prior=prior,
                        DElookup=DElookup)
# plotting and QC
# Principle Component Analysis on the counts
# plotPCA(Rcade)
# The MM plot shows log-ratios from DE plotted against log-ratios from the ChIP-seq. The colour of each
# point corresponds to the probability of both ChIP and DE being present
# plotMM(Rcade)
# The function plotBBB() plots log-odds values for ChIP-seq, DE, and combined ChIP-seq/DE analysis together.
# The package rgl is required for this 3D plot.
# plotBBB(Rcade)

# output
# Usually, the file of interest is “DEandChIP.csv", which contains the genes most likely to have both DE and
# ChIP signals. A full explanation of all of the files can be found in the help page: ?exportRcade
  out_dir <- sub(".csv", paste("/", as.character(j), sep = ""), output_dir[i])
  exportRcade(Rcade, directory=out_dir)
}
# By default, Rcade outputs the top 1000 geneIDs for each hypothesis. To increase the number of geneIDs,
# use a larger value for cutoffArg:
# exportRcade(Rcade, directory="RcadeOutput", cutoffArg=2000)

# Alternative cutoff methods are described in the exportRcade help page - for example,
# exportRcade(Rcade, directory="RcadeOutput", cutoffMode="B", cutoffArg=0)
}
