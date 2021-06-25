## General limma + edgeR analysis
library(tidyverse)
library(limma)
library(edgeR)


## Adding a function to make SAMPID SUBJID
## SAMPID to SUBJID function
SAMPID_to_SUBJID <- function(sampids){
    gregexpr("-" , "GTEX-1117F-0226-SM-5GZZ7")[[1]][[2]]
    
    greg_out <- gregexpr("-", sampids)
    lapply_out <- lapply(greg_out, function(x){x[[2]] -1})
    str_sub(sampids,1,lapply_out)
}

## Gene count location
gene_count_location <- "~/work2/tnieuwe1/data/gtex_v8/gtex-gene-counts-v8.rda"
## Subject data location
subject_location <-  "~/work2/tnieuwe1/data/gtex_v8/gtex_phenotypes_v8.csv"
## Matrisome data location
matrisome_location <- "~/work2/tnieuwe1/data/gtex_v8/matrisome_dat/matrisome_hs_masterlist_EDIT_r.csv"
## Age and sex output location
ans_output <- "~/work2/tnieuwe1/data/gtex_v8/matrisome_dat/age_sex_output/"
## Current tiss
cur_tiss <- "tish"

## Load in the read counts and other files
## read counts, gene names, sample data
load(gene_count_location)
## SUBJID data
subj_dat <- read.csv(subject_location)
## Matrisome data
mat_dat <- read.csv(matrisome_location)
colnames(mat_dat)[1] <- "Division" 

## Filter mat_dat
mat_dat_core <- filter(mat_dat, Division == "Core matrisome" | 
                           Category == "ECM Regulators")


##Filter to general and remove low count genes
ind <- stab$SMTSD %in% cur_tiss
stab_filt <- stab[ind,]
sub_dat <- dat[,ind]

## Filter out genes with low expression
keep.exprs <- filterByExpr(sub_dat)
sub_dat <- sub_dat[keep.exprs,]
gtab_filt <- gtab[keep.exprs,]


stab_filt$SUBJID <- SAMPID_to_SUBJID(stab_filt$SAMPID)
samp_n_subj <- left_join(stab_filt, subj_dat)
## Scale the age
samp_n_subj$SCALEAGE <- scale(samp_n_subj$AGE)
## Remove NA DTTHRDY
hrdy_ind <- !(is.na(samp_n_subj$DTHHRDY))
samp_n_subj_filt <- samp_n_subj[hrdy_ind,]
hrdy_ind_2 <- colnames(sub_dat) %in% samp_n_subj_filt$SAMPID 
sub_dat <- sub_dat[,hrdy_ind_2]

## Run Limma analysis
if (length(unique(samp_n_subj_filt$SEX)) > 1) {
    design <- model.matrix(~SCALEAGE + SEX +
                               ## adding new correction of batch effects
                               SMTSISCH +
                               as.factor(SMGEBTCH) +
                               factor(DTHHRDY,levels = c("1","2","3","4", "0"))
                           ,
                           data = samp_n_subj_filt)
} else{
    design <- model.matrix(~SCALEAGE +
                               ## adding new correction of batch effects
                               SMTSISCH +
                               as.factor(SMGEBTCH) +
                               factor(DTHHRDY,levels = c("1","2","3","4", "0"))
                           ,
                           data = samp_n_subj_filt)
}

## Normalize with voom
voom_dat <- voom(sub_dat, design = design)
## Run lmFit to fit the model
vfit <- lmFit(voom_dat, design)
## Run eBayes to generate statistics and DE results.
efit <- eBayes(vfit)

### Clean up data
## Pull out the results of SCALEDAGE
age_results_limma <- topTable(efit,coef=2, sort = "none", n =Inf)
age_results_mat <- age_results_limma %>%
    rownames_to_column(var = "gene_id") %>%
    left_join(., gtab_filt) %>%
    filter(gene_name %in% mat_dat_core$Gene.Symbol)

## Pull out the results of SEX
if (length(unique(samp_n_subj$SEX)) > 1) {
    sex_results_limma <- topTable(efit,coef=3, sort = "none", n =Inf)
    sex_results_mat <- sex_results_limma %>%
        rownames_to_column(var = "gene_id") %>%
        left_join(., gtab_filt) %>%
        filter(gene_name %in% mat_dat_core$Gene.Symbol)
}



## save output
write.csv(age_results_mat, paste0(ans_output, "limma_general_res_age_batch.csv"),row.names = F)
if (length(unique(samp_n_subj$SEX)) > 1) {
    write.csv(sex_results_mat, paste0(ans_output, "limma_general_res_sex_batch.csv"),row.names = F)
}