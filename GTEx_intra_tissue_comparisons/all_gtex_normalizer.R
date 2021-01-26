#Generate a variance stabilized transformed dataframe

source("~/work2/tnieuwe1/data/gtex_v7/tim_scripts.R")
library(DESeq2)
library(tidyverse)


###load datasets (heavy part of code)
#load("~/data2/gtex_v8/gtex-gene-counts-v8.rda")
load("~/work2/tnieuwe1/data/gtex_v8/gtex-gene-counts-v8.rda")
matrisome <- read.csv("~/work2/tnieuwe1/data/gtex_v8/matrisome_dat/matrisome_hs_masterlist_EDIT_r.csv")

generaltab <- stab
generaldat <- dat

## load subject annotation
# subjtab <- read.table("~/data2/gtex_v8/GTEx_v8_Annotations_SubjectPhenotypesDS.txt",
#                       header=TRUE,stringsAsFactors=FALSE,sep="\t",fill=TRUE,
#                       quote="")

subjtab <- read.table("~/work2/tnieuwe1/data/gtex_v8/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
                      header=TRUE,stringsAsFactors=FALSE,sep="\t",fill=TRUE,
                      quote="")

##Create shortnames to the data and attatch phenotype data

phen <- read.csv("~/work2/tnieuwe1/data/gtex_v8/matrisome_dat/phenotypes_light_gtex_v8.csv", stringsAsFactors = F)
colnames(phen)[1] <- "SAMPID"

#Below merge on subject ID
finaltab <- left_join(generaltab, phen, by = "SAMPID")

rownames(generaltab) <- generaltab$SAMPID

#Create the Deseq 2 object
generalDES <- DESeqDataSetFromMatrix(countData = generaldat,
                                     colData = finaltab,
                                     rowData = gtab,
                                     design = ~ 1)
rownames(generalDES) <- gtab$Description


## variance stabilizing transformation
matrisomeVSD <- vst(generalDES, blind = T, nsub = 1000, fitType = "parametric")

#Create matrisome index
mat_ind <- gtab$Description %in% matrisome$Gene.Symbol

#Save the normalized output
save(gtab, matrisomeVSD, finaltab, mat_ind, file="~/work2/tnieuwe1/projects/matrisome/nmatrisomeVSD.rda")
