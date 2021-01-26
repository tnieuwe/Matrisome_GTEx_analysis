library(DESeq2)
library(gplots)
library("pheatmap")
library("RColorBrewer")
library(SuppDists)
library(gplots)
library(fields)
library(GMD)
library(tidyverse)



###general Anaylsis


library(DESeq2)
#setwd("~/Dropbox/projects/gtex/artery")
load(file="~/work2/tnieuwe1/data/gtex_v8/diff_samps/general/general-vsd-mean-filtered.rda")

###MAKE general GENERAL FOR TEST
VSDMeanFiltered <- generalVSDMeanFiltered

##SUBSET OUT MATRISOME
matrisome <- read.csv("~/work2/tnieuwe1/projects/matrisome/matrisome_hs_masterlist_EDIT_r.csv", stringsAsFactors = F)
colnames(matrisome)[1] <- "division"

matrisome_core <- filter(matrisome, division == "Core matrisome" | Category == "ECM Regulators")
matrisome_assoc <- filter(matrisome, division == "Matrisome-associated")

ind <- (gtabMeanFiltered$gene_name %in% matrisome_core$Gene.Symbol)

gtab_mat <- gtabMeanFiltered[ind,]

#Get matrix out of data
vsdMeanFiltered <- generalVSDMeanFiltered

mat_dat <- assay(vsdMeanFiltered)

mat_dat <- mat_dat[ind,]

####

vsdMeanFiltered <- generalVSDMeanFiltered

## calculate variance for each genes
v <- apply(mat_dat,1,var)
ind <- which(v > 1.5)
names(v) <- gtab_mat$gene_name
v <- as.data.frame(v)
v <- rownames_to_column(v, var = "Gene.Symbol")
v <- v %>% mutate(v3 = v > 3, v2 = v >2) %>% mutate(variance = if_else(v2 == T, if_else(v3 == T, "V3", "V2"), "V0" ))



## heatmap of high variance genes
library(gplots)
#creating highvar subset
vsdHighVar <- mat_dat[ind,]
rownames(vsdHighVar) <- gtab_mat$gene_name[ind]
colnames(vsdHighVar) <- vsdMeanFiltered$SAMPID
vsdHighVarCentered <- vsdHighVar - rowMeans(vsdHighVar)
pdf("~/work2/tnieuwe1/data/gtex_v8/matrisome_dat/mat_cluster_output/general-matrisome-high-variance-heatmap.pdf", width=8, height=8)
hm <- heatmap.2(vsdHighVarCentered, trace="none",col=bluered(20), 
                scale="row", margins=c(4,4), cexCol=0.25, cexRow=0.25,
                density.info="histogram",key.title="")
#dev.off()

## correlation between high variance genes
library("pheatmap")
library("RColorBrewer")

#Removing duplicated genes, Tim code
dim(vsdHighVar)
vsdHighVar <-  vsdHighVar[!(duplicated(rownames(vsdHighVar))),]
dim(vsdHighVar)

#Gcor is correlations between each genes
gcor <- cor(t(vsdHighVar), method="kendall")

## critical value
#library(pvrank)
library(SuppDists)
ngene <- nrow(vsdHighVar)
nsamp <- ncol(vsdHighVar)
#tauCrit_old <- round(-qrank(0.01/((ngene^2)-ngene)/2, nsamp, index="kendall", approx = "gaussian")$Cq, 2)
tauCrit <- round(qKendall(0.01/((ngene^2)-ngene)/2, nsamp, lower.tail = F), 2)

## cluster based on absolute correlation distance
dendo <- hclust(d=as.dist(1-(gcor)), method = "average")

## clean up heatmap to focus on interesting gene clusters
clusts <- cutree(dendo, h=1-tauCrit)
tabclusts <- table(clusts)
#arbitrary piece of code that allows usd to select miminum gene in a cluster
ind <- which(tabclusts >= 6) #arbitrary
clusts[!(clusts %in% names(tabclusts)[ind])] <- 0
ntmp <- names(clusts)
clusts <- as.numeric(as.factor(clusts))-1
names(clusts) <- ntmp

## heatmap

#Naming clusters based on their clust values
#Very manual process
n <- length(table(clusts))
cnames <- rep("", length(clusts))

#Tim for loop
for (clst in 1:(n)){
  cnames[clusts==clst] <- LETTERS[clst]
}

#Addition of matrisome or matrisome associated code here
anno = data.frame("cluster"=factor(cnames))
rownames(anno) = names(clusts)

#Anno 2
#gene_4_anno <- rownames(anno)
#gene_4_anno <- as.data.frame(gene_4_anno)
#colnames(gene_4_anno) <- "Gene.Symbol"
#anno_2 <- left_join(gene_4_anno, matrisome) %>% select(division)
#rownames(anno_2) <- gene_4_anno$Gene.Symbol

#Anno 2 variance
gene_4_anno <- rownames(anno)
gene_4_anno <- as.data.frame(gene_4_anno)
colnames(gene_4_anno) <- "Gene.Symbol"
anno_2 <- left_join(gene_4_anno, v) %>% select(variance)
rownames(anno_2) <- gene_4_anno$Gene.Symbol


#anno_colors = list("cluster"=c("white", brewer.pal(n-1, "Spectral")))
anno_colors = list("cluster"=c("white", brewer.pal(n-1, "Spectral")))
names(anno_colors$cluster) <- c("", LETTERS[1:(n)-1])
colors <- colorRampPalette( brewer.pal(11, "RdBu")[2:10] )(255)
brks <- seq(0,2,length=length(colors)+1)

chm <- pheatmap(1-gcor, col=colors, breaks=brks,
                cluster_rows=dendo, cluster_col=dendo, 
                legend_breaks=c(2,1.5,1,0.5,0), 
                legend_labels=c("-1.0","-0.5","0","0.5","1.0"),
                main="", fontsize=5, fontsize_row=5, fontsize_col=5,
                
                annotation_col = anno, annotation_row = anno_2, 
                
                annotation_colors = anno_colors,
                border_color = NA,
                filename="~/work2/tnieuwe1/data/gtex_v8/matrisome_dat/mat_cluster_output/general-between-gene-correlation-high-variance-genes_mat.pdf", 
                width=10, height=8)
#dev.off()

## make matrix to output
## make matrix to output
gclusts <- matrix("", nrow=max(table(cnames)), ncol=length(table(cnames)))
colnames(gclusts) <- levels(anno$cluster)
for(k in 1:length(table(cnames))){
  tmp <- colnames(gclusts)[k]
  gclusts[1:sum(cnames==tmp), k] <- rownames(gcor)[cnames==tmp]
}
#Generate sub clusters
for (column in LETTERS[1:(length(unique(cnames)) -1 )]) {
  #Pull out the first column
  test_col <- gclusts[,column]
  first_gene <- test_col[1]
  #Make an index to subset gcor to just a single clusters
  ind <- which(rownames(gcor) %in% test_col)
  sub_cor <- gcor[ind,ind]
  #Test if the cluster has any negatively correlating genes
  ind2 <- sub_cor[,first_gene] < 0
  
  #If there are negatively correlating genes...
  if (sum(ind2) > 0) {
    #generate sub cluster 1
    sub_1 <- rownames(sub_cor)[!ind2]
    #Subcluster 2
    sub_2 <- rownames(sub_cor)[ind2]
    
    #Remover old column to replace with new columns
    gclusts <- as.matrix(gclusts[,!(colnames(gclusts) == column)])
    
    
    #nrow(gclusts) Not sure if useful, keeping for now
    #add sub-columns
    
    #Make a single matrix of same rows as gclusts, one for each sub-cluster. 
    matrix_1  <- as.matrix(c(sub_1 ,rep("", nrow(gclusts) - length(sub_1))), )
    colnames(matrix_1) <- paste0(column,"-1")
    
    # cbind(gclusts, matrix_1, matrix_2)
    
    
    matrix_2  <- as.matrix(c(sub_2 ,rep("", nrow(gclusts) - length(sub_2))), )
    colnames(matrix_2) <- paste0(column,"-2")
    
    #Add subclusters to gclusts
    gclusts <-   cbind(gclusts, matrix_1, matrix_2)
    
  }
  
  
}



#Realphebatize columns for sake of beauty 
gclusts <- gclusts[,order(colnames(gclusts))]
write.csv(gclusts, file="~/work2/tnieuwe1/data/gtex_v8/matrisome_dat/mat_cluster_output/general-gene-clusters-high-variance_mat.csv",
          quote=FALSE, row.names=FALSE)

##############

## stratify data by inclusion in clusters 
ind <- which(rownames(vsdHighVar) %in% rownames(gcor)[clusts>0])
vsdHighVarClust <- vsdHighVar[ind,]
vsdHighVarCenteredClust <- vsdHighVarCentered[ind,]
vsdHighVarNoise <- vsdHighVar[-ind,]

gcorClust <- gcor[ind,ind]
dendoClust <- hclust(d=as.dist(1-(gcorClust)), method = "average")
cnamesClust <- cnames[ind]

## clusters A-H
n <- length(table(cnamesClust))
anno = data.frame("cluster"=factor(cnamesClust))
rownames(anno) = rownames(vsdHighVarClust)
#Getting annotations for matrisome from previous samples
anno_2_clean <- subset(anno_2, rownames(anno_2) %in% rownames(anno))
anno_colors = list("cluster"=brewer.pal(n, "Spectral"))
names(anno_colors$cluster) <- LETTERS[1:n]
colors <- colorRampPalette( brewer.pal(11, "RdBu")[2:10] )(255)
brks <- seq(0,2,length=length(colors)+1)
chm <- pheatmap(1-gcorClust, col=colors, breaks=brks,
                cluster_rows=dendoClust, cluster_col=dendoClust, 
                legend_breaks=c(2,1.5,1,0.5,0), 
                legend_labels=c("-1.0","-0.5","0","0.5","1.0"),
                main="", fontsize=5, fontsize_row=3, fontsize_col=3,
                annotation_col = anno, annotation_row = anno_2_clean, 
                annotation_colors = anno_colors,
                border_color = NA,
                filename="~/work2/tnieuwe1/data/gtex_v8/matrisome_dat/mat_cluster_output/general-between-gene-correlation-high-variance-genes-clean-mat.pdf", 
                width=10, height=8)
#dev.off()

## remake with non A-H genes
gcorNoise <- gcor[-ind,-ind]
dendoNoise <- hclust(d=as.dist(1-(gcorNoise)), method = "average")
geneNoise <- dendoNoise$labels
anno_2_noise <- subset(anno_2, rownames(anno_2) %in% geneNoise)
colors <- colorRampPalette( brewer.pal(11, "RdBu")[2:10] )(255)
brks <- seq(0,2,length=length(colors)+1)
chm <- pheatmap(1-gcorNoise, col=colors, breaks=brks,
                cluster_rows=dendoNoise, cluster_col=dendoNoise, 
                legend_breaks=c(2,1.5,1,0.5,0), 
                legend_labels=c("-1.0","-0.5","0","0.5","1.0"),
                main="", fontsize=10, fontsize_row=10, fontsize_col=10,
                annotation_row = anno_2_noise, 
                border_color = NA,
                filename="~/work2/tnieuwe1/data/gtex_v8/matrisome_dat/mat_cluster_output/general-between-gene-correlation-high-variance-genes-noise-mat.pdf", 
                width=10, height=8)
#dev.off()


#####################

## add subject specific data
subjtab <- read.table(file="~/work2/tnieuwe1/data/gtex_v8/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", fill=T,
                      stringsAsFactors=FALSE, sep="\t", header=T)
subjIDs <- sapply(strsplit(colnames(vsdHighVarClust),"-"),function(x) paste(x[1],x[2],sep="-"))
map <- match(subjIDs,subjtab$SUBJID)
identical(subjIDs,subjtab$SUBJID[map])

gender <- subjtab$SEX[map]
age <- subjtab$AGE[map]
death <- subjtab$DTHHRDY[map]
death3 <- death
death3[death3>1] <- 2

## add sample specific data
samptab <- read.table("~/work2/tnieuwe1/data/gtex_v8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
                      header=TRUE,stringsAsFactors=FALSE,sep="\t",fill=TRUE,
                      quote="")
map <- match(colnames(vsdHighVarClust),samptab$SAMPID)
identical(colnames(vsdHighVarClust),samptab$SAMPID[map])
samptab <- samptab[map,] 

## extract batch date as a numeric factor
batch <- as.Date(samptab$SMNABTCHD, format="%m/%d/%Y")
batch2 <- as.numeric(as.factor(batch))
bcols <- colorRampPalette( brewer.pal(9, "Purples")[2:9] )(90)

## extract RNA integrity number
rin <- samptab$SMRIN
rcols <- colorRampPalette( brewer.pal(9, "YlGn"))(30)

#Row colors
#make anno_colors data frame and add to annotatiion
#map3_color
map3_anno <- anno
map3_anno <- rownames_to_column(anno, var = "gene")

#Remove NA if there
if (is.na(names(anno_colors$cluster[length(anno_colors$cluster)])) ) {
  anno_colors  <- anno_colors$cluster[-length(anno_colors$cluster)]
}


map3_color <- as.data.frame(anno_colors)
colnames(map3_color) <- "color"
map3_color <- rownames_to_column(map3_color, var = "cluster")
#Make matrisome labels for map3

#Generating colors to be used and data frame required to be added into rlab
colors <- as.data.frame(cbind(c("Matrisome-associated","Core matrisome"), (RColorBrewer::brewer.pal(3, "Paired")[-3])))
colnames(colors) <- c("division", "Division")
anno_mat <- left_join(gene_4_anno, matrisome) %>%
  left_join( colors)
anno_mat <- anno_mat %>% select(Gene.Symbol, Division)
colnames(anno_mat)[1] <- "gene"

map3_full <- left_join(map3_anno, map3_color) %>% left_join(anno_mat)
map3_full <- map3_full[order(map3_full$cluster, decreasing = T),]


#Making age a thing with r color brewer

##No bin
#age_color <- brewer.pal(6, "Oranges")
#age_df <- as.data.frame(cbind(c("20-29", "30-39", "40-49", "50-59", "60-69","70-79"), age_color ))

##With bin
age_color <- brewer.pal(5, "Oranges")
age_df <- as.data.frame(cbind(c("20-39", "40-49", "50-59", "60-69","70-79"), age_color ))

colnames(age_df) <- c("age","color")
age <- as.data.frame(age, stringsAsFactors = FALSE)
#bin
age$age[age$age == "20-29"] <- "20-39"
age$age[age$age == "30-39"] <- "20-39"
age_full <- left_join(age, age_df)


## heatmap with sample annotation
library(gplots)
library(fields)
library(GMD)
library("devtools")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

#source("heatmap.3.R")
clab <- cbind(rcols[cut(rin, 30)],
              bcols[batch2], 
              c("grey","black")[gender],
              c("dark green","green","light green")[death3+1],
              as.character(age_full$color)
              
)
colnames(clab) <- c("RIN","Date","Gender","Death", "Age")

map3_full <- select(map3_full, color, Division) 
#rlab <- t(as.matrix(map3_full$color))
rlab <- t(as.matrix(map3_full))

pdf("~/work2/tnieuwe1/data/gtex_v8/matrisome_dat/mat_cluster_output/general-mat-high-var-heatmap-clean-subj-info.pdf", width=9.5, height=8)
hm3 <- heatmap.3(vsdHighVarCenteredClust[order(cnamesClust, decreasing = TRUE),], 
                 Rowv=NA, trace="none",
                 col=bluered(20),
                 scale="none", 
                 margins=c(2,12), labCol=NA, labRow=NA,
                 density.info="none", KeyValueName="Residuals", keysize=1,
                 xlab="tish Samples", ylab="", dendrogram="none",
                 sidesize=4, 
                 ColSideColors=clab, ColSideColorsSize=4,
                 RowSideColors = rlab
                 #add.exprs=mtext(side=2,paste("Cluster", LETTERS[1:8]), at=c(0.15,0.47,0.68,0.8,0.88,0.95,1.01),line=2, las=2)
)
##dev.off()

legend(x=0.86,y=0.97,c("Ventilator","Violent","Natural"), title="Death",
       fill=c("dark green","green","light green"), border=FALSE, bty="n", cex=1)
legend(x=0.86,y=0.81,c("Female","Male"), title="Gender",
       fill=c("black","grey"), border=FALSE, bty="n", cex=1)
legend(x=0.86,y=0.68,c("06/27/11","","09/12/13","","01/28/14"), title="Date",
       fill=rep("white",5), border=FALSE, bty="n", cex=1)
colorbar.plot(x=0.89, y=0.55, strip=sort(unique(batch2), decreasing=TRUE), col=bcols, 
              strip.width=0.02, strip.length=0.14, horizontal=FALSE)
legend(x=0.86,y=0.45,c("3.2","","7.0","","9.4"), title="  RIN",
       fill=rep("white",5), border=FALSE, bty="n", cex=1)
colorbar.plot(x=0.89, y=0.32, strip=seq(from=max(rin), to=min(rin), length=length(rin)), 
              col=rcols, strip.width=0.02, strip.length=0.14, horizontal=FALSE)
#Tim added cluster bar

legend(x=0.86,y=0.22,c("20-39", "40-49", "50-59", "60-69","70-79"), title="Age",
       fill=age_color, border=FALSE, bty="n", cex=1)

legend(x=0.01,y=0.5,map3_color$cluster, title="Cluster",
       fill=as.character(map3_color$color), border=FALSE, bty="n", cex=1)

legend(x=0.01,y=0.3,c("Core","Assoc"), title="Matrisome",
       fill=c("#1F78B4","#A6CEE3"), border=FALSE, bty="n", cex=1)
#dev.off()


save(vsdHighVar, vsdHighVarCentered, vsdHighVarClust, vsdHighVarNoise, 
     vsdHighVarCenteredClust, cnames, cnamesClust,
     file="~/work2/tnieuwe1/data/gtex_v8/diff_samps/general/general-high-variance-gene-data-mat.rda")

#########################################

## compute avg cluster expression for each sample
load(file="~/work2/tnieuwe1/data/gtex_v8/diff_samps/general/general-high-variance-gene-data-mat.rda")

## switch direction of negatively correlated gene
#If statement added by Tim for sex based collections
if ("XIST" %in% vsdHighVarCenteredClust) {
  vsdHighVarCenteredClust["XIST",] <- -vsdHighVarCenteredClust["XIST",]
}

## extract expression for each cluster, convert to zscores, and summarize profile
clusterProfiles <- matrix(nrow=ncol(vsdHighVarCenteredClust), ncol=length(unique(cnamesClust)))
for(k in 1:ncol(clusterProfiles)){
  ind <- which(cnamesClust == LETTERS[k])
  etmp <- vsdHighVarCenteredClust[ind,]
  ztmp <- etmp / apply(etmp, 1, mad)
  clusterProfiles[,k] <- colMeans(ztmp)
}
rownames(clusterProfiles) <- colnames(vsdHighVarCenteredClust)
colnames(clusterProfiles) <- LETTERS[1:ncol(clusterProfiles)]

save(clusterProfiles, file="~/work2/tnieuwe1/data/gtex_v8/diff_samps/general/general-cluster-profiles-mat.rda")
write.csv(clusterProfiles, file="~/work2/tnieuwe1/data/gtex_v8/matrisome_dat/mat_cluster_output/general-cluster-profiles.csv", quote=FALSE)
write.csv(v, file="~/work2/tnieuwe1/data/gtex_v8/matrisome_dat/mat_cluster_output/general-matrisome-variance.csv", quote=F, row.names = F)
