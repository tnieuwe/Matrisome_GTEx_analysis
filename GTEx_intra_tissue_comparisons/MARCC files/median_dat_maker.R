## This file was ran on the MARCC cluster

library(DESeq2)
library(tidyverse)
library(matrixStats)
load(nmatrisomeVSD.rda)
dat <- assay(matrisomeVSD)
tiss_list <- unique(finaltab$SMTSD)

median_dat <- NULL
i <- 1
for (tiss in tiss_list) {
  ind <- finaltab$SMTSD %in% tiss
  sub_dat <- dat[,ind]
  new_median <- rowMedians(sub_dat)
  median_dat <- cbind(median_dat, new_median)
  colnames(median_dat)[i] <- tiss
  
  i <- i + 1
}
save(median_dat, gtab, file ="median_gtex_v8.rda")