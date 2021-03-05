#Tim script


###Options
options(stringsAsFactors = FALSE)





#####General

#Name: Look()

#Purpose: To look at the top left corner of the data frame, somewhat like head, but with a colname() limit.

look <- function(dat, row_num = 0, col_num = 0){
  
  if(row_num == 0 || row_num > dim(dat)[1]){
    row_num <- ifelse(dim(dat)[1] < 5, dim(dat)[1], 5)
  }
  
  if(col_num == 0 ||  col_num > dim(dat)[2] ){
    
    col_num <- ifelse(dim(dat)[2] < 5, dim(dat)[2], 5)
    
  }
  
  dat[1:row_num, 1:col_num]
}



####GTEx Manipulation

#Name: shorter_gtex_names_df()

#Purpose: To change long Gtex names into the names of the individuals the tissue came from. This version specifically returns a dataframe

shorter_gtex_names_df <- function(x, col =1 ){
  col2 <- colnames(x)[col]
  #If else statement used to determine if the name is long or short, if already short nothing should happen, if long, then it will be appropriately shortened. This code works off the assumption that the gtex IDs are either 4 or 5 characters in length
  if(nchar(x[1, col]) > 10){
    
    #Generating stand in list
    shortname <- vector()
    #Used subset to create a vector, drop needed to be true to make the tibble output a vector
    for (k in subset(x, select =1, drop = T)) {
      #If else for 5 or 4 length IDs
      #Using the "-" to determine if the id has 4 or 5 characters
      if(substr(k, start = 10, stop = 10) == "-"){
        #4 length
        m <- substr(k, 1, 9)
        shortname <- c(shortname, m)
      } else{
        #5 length
        m <- substr(k, 1, 10)
        
        shortname <- c(shortname, m)
        
      }
      
    }
    x[col] <- shortname
    
  } else{
    #Else just puts the shortnames back in case they are already short, thus making the code safe to run on all text.
    shortname <- x[col]
    x[col] <- shortname
  }
  return(x)
}

#Name: shorter_gtex_names_vect()

#Purpose: To change long Gtex names into the names of the individuals the tissue came from.

shorter_gtex_names_vect <- function(x, col =1 ){
  col2 <- colnames(x)[col]
  #If else statement used to determine if the name is long or short, if already short nothing should happen, if long, then it will be appropriately shortened. This code works off the assumption that the gtex IDs are either 4 or 5 characters in length
  if(nchar(x[1, col]) > 10){
    
    #Generating stand in list
    shortname <- vector()
    #Used subset to create a vector, drop needed to be true to make the tibble output a vector
    for (k in subset(x, select =1, drop = T)) {
      #If else for 5 or 4 length IDs
      #Using the "-" to determine if the id has 4 or 5 characters
      if(substr(k, start = 10, stop = 10) == "-"){
        #4 length
        m <- substr(k, 1, 9)
        shortname <- c(shortname, m)
      } else{
        #5 length
        m <- substr(k, 1, 10)
        
        shortname <- c(shortname, m)
        
      }
      
    }
    
    
  } else{
    #Else just puts the shortnames back in case they are already short, thus making the code safe to run on all text.
    shortname <- x[col]
  
  }
  return(shortname)
}

#Name: age_to_num

#Purpose: Make ages standardized numbers

age_to_num <- function(x){
  
  k <- x$AGE
  x$AGE <- as.integer(substr(k, start = 1, stop = 2))+5
  return(x)
}



####Gtex tests

#Name: correlation_graph()

#Purpose: Give correlation output along with a graph of the data
correlation_graph <- function(dats, in_var, de_var, cortest = T, graph = "point", corre = "pearson", lm = F){

  if(lm == F){
    plot <- ggplot(data = dats, aes_string(x = in_var , y = de_var, group = in_var))
  }
  if(lm == T){
    plot <- ggplot(data = dats, aes_string(x = in_var , y = de_var))
  }
  if(graph == "violin"){
    grph <- plot + geom_violin()
  }
  if(graph == "boxplot"){
    grph <-  plot + geom_boxplot()
  }
  if(graph == "point"){
    grph <-  plot + geom_point() 
  }
  
  if(cortest == T){
    test.res <- cor.test(dats[[in_var]], dats[[de_var]] , method = corre)
    return(list(grph, test.res)) 
  } else {
    return(grph)
  }
  
}

#Name: correlation_all_sig
#Purpose: Finding a correlation between a Z score and all genes.

correlation_all_sig <- function(dats, fixed, start_col = 4, end_col = ncol(dats), meth = "pearson", sig_only = T, pval_meth = "bonferroni"){
  #Run correlations, pull out pval
  out1 <- sapply(dats[,start_col:end_col], function(i) cor.test(dats[[fixed]], i, method = meth)$p.value)
  out_adjust <- p.adjust(out1, method = pval_meth)
  sig <- na.omit(out_adjust[out_adjust<0.05])
  sig_nocorrect <- na.omit(out1[out_adjust<0.05])
  #Pull sig names for cluster analysis
  sig_names <- names(sort(sig_nocorrect))
  
  #Correlation estimate
  out2 <- sapply(dats[,start_col:end_col], function(i) cor.test(dats[[fixed]], i, method =meth)$estimate)
  
  #String_i removing the endlines
  names(out2) <- str_sub(names(out2), 1, -5)        
  
  sig2 <- out2[names(out2) %in% names(sig_nocorrect)]
  
  if(sig_only == T){
    stats_output <- data.frame(sig_nocorrect, sig, sig2)
    stats_sort <- stats_output[order(sig_nocorrect),]
    colnames(stats_sort) <- c("p-value", "correct p-val", "correlation")
    return(stats_sort)
  }      else{
    stats_output <- data.frame(out1, out_adjust, out2)
    stats_sort <- stats_output[order(out1),]
    colnames(stats_sort) <- c("p-value", "correct p-val", "correlation")
    return(stats_sort)
  }
}


#Name: t_test_graph()
#Purpose: Gives a t-test output along with a graphical representation of a t-test


z_test_graph <- function(dats, in_var, de_var, test = "norm", graph = "vio_jitter"){
  
  library(BSDA)
  
  dats <- dats %>% filter( UQ(as.name(in_var)) != 99)
  
  plot <- ggplot(data = dats, aes_string(x = in_var , y = de_var, group = in_var))
  
  
  if(graph == "vio_jitter"){
    grph <- plot + geom_violin() + geom_jitter(width = 0.4)
  }
  if(graph == "violin"){
    grph <- plot + geom_violin()
  }
  if(graph == "boxplot"){
    grph <-  plot + geom_boxplot()
  }
  if(graph == "point"){
    grph <-  plot + geom_point()
  }
  if(graph == "jitter"){
    grph <-  plot + geom_jitter()
  }
  
  if(test == "norm"){
    test.res <- z.test(as.formula(paste(de_var, in_var, sep = "~")), data = dats)
    return(list(grph, test.res)) 
  } 
  if(test == "pair"){
    test.res <- pairwise.t.test(dats[[de_var]],dats[[in_var]],p.adjust="bonf")
    return(list(grph, test.res)) 
  } else {
    return(grph)
  }
  
}

#Name: t_test_all_sig
#Purpose: Finding a correlation between a Z score and all genes.

z_test_all_sig <- function(dats, fixed, start_col = 2, end_col = ncol(dats), sig_only = T, pval_meth = "hochberg"){
  library(BSDA)
  #Note the use of UQ, it is a dplyr script that allows the "quoted" column to be used as a filterable function
  dats <- dats %>% na.omit() %>% filter( UQ(as.name(fixed)) != 99)
  
  
  #Below code removes columns that all have the same values
  if(!(length(colnames(dats)) == length(colnames(dats[vapply(dats, function(x) length(unique(x)) > 1, logical(1L))])))){
    dif <- length(colnames(dats)) - length(colnames(dats[vapply(dats, function(x) length(unique(x)) > 1, logical(1L))]))
    end_col <- end_col - 1
    dats <- dats[vapply(dats, function(x) length(unique(x)) > 1, logical(1L))]
    
  }
  
  #Set as factor
  dats[[fixed]] <- as.factor(dats[[fixed]])
  
  out <- sapply(dats[,start_col:end_col], function(i) z.test(i ~ dats[[fixed]])$p.value)
  
  out_adjust <- p.adjust(out, method = pval_meth)
  
  
  out2 <- sapply(dats[,start_col:end_col], function(i) z.test(i ~ dats[[fixed]])$estimate)
  out2 <- t(out2)
  out2 <- data.frame(out2)
  
  
  sig <- out_adjust[out_adjust<0.05]
  sig_nocorrect <- out[out_adjust<0.05]
  
  sig2 <- out2[rownames(out2) %in% names(sig_nocorrect),]
  
  sig_names <- names(sort(sig_nocorrect))
  
  if(sig_only == T){
    stats_output <- data.frame(sig_nocorrect, sig, sig2)
    stats_sort <- stats_output[order(sig_nocorrect),]
    colnames(stats_sort) <- c("p-value", "correct p-val","mean group 0", "mean group 1")
    return(stats_sort)
  }      else{
    stats_output <- data.frame(out, out_adjust, out2)
    stats_sort <- stats_output[order(out),]
    colnames(stats_sort) <- c("p-value", "correct p-val","mean group 0", "mean group 1")
    return(stats_sort)
  }
}


#Name: in_clust_df()
#Purpose: Find if significant genes are inside dataframes, specificaly used on cluster output files

#Function testing for shared genes in dataframes

in_clust_df <- function(sig, clusts){
  
  in_cluster <- vector()
  clusts <- unlist(clusts)
  sig_1 <- sig[,1]
  for(i in sig_1){
    if(i %in% clusts){
      in_cluster <- c(in_cluster, names(clusts)[i == clusts] )
    }else{
      in_cluster  <- c(in_cluster, NA)
    }
  }
  data.frame(sig, in_cluster) 
  
}

#Name: usable_gtex
#Purpose: Turns gtex data output into a more manageable format to be used in analysis


usable_gtex <- function(dat, tab){
  #Exract matrix as dataframe
  matrx <- as.data.frame(assay(dat))
  #Pull ens gene names so they can be replaced with gene-ids from the data table
  matrx$ens <- row.names(matrx)
  matrx_filt <- matrx %>% filter(ens %in% gtabMeanFiltered$Name) %>% mutate(
    gene_id = gtabMeanFiltered$Description) %>% select(
      ens, gene_id, everything()
    )
  #Remove ens and transform for analysis
  matrx_filt <-   matrx_filt[,-1]
  gene_ids <-  matrx_filt$gene_id
  t_genes_dat <- t( matrx_filt[-1])
  colnames(t_genes_dat) <- gene_ids
  dat2 <- as.tibble(t_genes_dat, rownames = "gtex-id")
  
  
  #Shorten name length
  dat_fin <- shorter_gtex_names(dat2)
  return(dat_fin)
}
