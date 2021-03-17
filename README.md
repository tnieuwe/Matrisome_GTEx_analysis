# Matrisome_GTEx_analysis
The Github Repository for the paper "Sex, age, tissue, and disease patterns of matrisome expression in GTEx transcriptome data" by Tim O. Nieuwenhuis, Avi Z. Rosenberg, Matthew N. McCall, and Marc K. Halushka. 

## Purpose
This Github Repo exists for the purpose of giving individuals the ability to reproduce the findings and figures found in the previously named paper. All code was written by Tim Nieuwenhuis and as a result any issues with the Github should be directed towards him via issues or email correspodence. 

## Organization
* **global_in**: Holds files used in multiple analyses
    *  matrisome_hs_masterlist_EDIT_r.csv: A csv version of the list of Matrisome Genes made by the Hynes lab. Original data can be retrieved [here](http://matrisomeproject.mit.edu/static/media/uploads/Logos/.thumbnails/excel-2010-icon.jpg/excel-2010-icon-40x30.jpg).
    *  **tim_scripts.R**: A file holding small quality of life scripts that are used when working with GTEx data.
    *  **gtex_phenotypes_v8.csv** *not available*: This is a file that includes specific phenotype data related to the individuals in GTEx such as their true age. This data is not included in the git hub as it is only accessible through DBGap with the correct permissions. The main usage of this tool is for non-binned ages of the individuals in GTEx, however one can use the binned ages available [here](https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.xlsx) as a substitute. 
    *  **GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt**: Publically available sample data from GTEx
    *  **GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt**: Publically available subject data from GTEx
* **GTEx_protein_gene_correlations**: Looking at GTEx Specific protein gene correlations and correlations from Wang et al. Includes code for Figure 1A, Figure 1B, Figure S1, and Figure S2.
  * **protein_gene_correlation_analysis.rmd**: This file includes the analysis of correlations from [Jiang et al.](https://www.sciencedirect.com/science/article/pii/S0092867420310783). The associated data files from the paper required to run the code are to large to upload to Github, and so their https are found in the comments. This file includes code for Figure S1.
  * **wang_et_al_correlations.Rmd**: This file includes the analysis on [Wang et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6379049/)'s paired protein gene analysis of tissues. This file is used to generate Figure 1A, 1B, and Figure S2. 
* **GTEx_intra_tissue_comparisons**: Files and data used to generate Figure 1C and. 1D
  * **MARCC files**: The files used to generate the median GTEx data used in this analysis. Note these files were not ran locally but on the MARCC cluster.
      * **all_gtex_normalizer.R**: The code used to normalize all of the GTEx data together
      * **median_dat_maker.R**: The code used to  to generated a median value for reach tissue from the previous code.
   * **median_matrisome_analysis_regulator_updated.rmd**: Analyzed the median GTEx data. This analysis looks into how different tissues cluster together based on genes found in the matrisome and what percentage of their expression is from the matrisome. Generates Figure 1C and 1D.
   * **input**:
      * gtex_name_translator: Used to connect GTEx tissue names through their various forms in this dataset (normal string, snake_case, abbreviations).
      * median_gtex_v8: The output from the MARCC files that includes all tissues and their median VST normalized expression.
      * median_tpm_v8: The median TPM data from the GTEx portal, download [here](https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz)
* **GTEx_TCGA_comparison**: The analysis comparing GTEx normal tissue and TCGA normal tissue, against TCGA cancer tissue. Includes code for Fig 4.
   * **quantile_rank_analysis**: An rmd that allows one to go through the entire analysis for whichever tissue they are interested in.
   * **quantile_rank_analysis_all_tissue**: An R file that is a loop of **quantile_rank_analysis** going through all tissues and generating all required files for **rank_output_analysis**.
   * **rank_output_analysis**: The analysis that compares the various rank changes between normal (GTEx and TCGA) and cancer (TCGA) tissue. Generates Fig 4.
* **GTEx_IPF_analysis**: An analysis using an IPF experiment and normal GTEx lung data to see if GTEx can be used as normal tissue in this context. Includes code for Fig 5A, 5B, Table 1, and Fig S5.
   * **GSE134692_ipf_cluster_maker**: A file used to generate high variance clusters in the IPF dataset GSE134692.
   * **GSE134692_ipf_analysis_w_gtex**: A file, using clusters generated in the previous file to test if GTEx's samples will cluster with normal tissues from the high variance genes associated with IPF.
   * **data_in**:
      * **gtex_pneumo**: Data classifying GTEx lung sample state as defined by a pathologist.
      * **GSE134692**: Alld ata from IPF experiment available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134692
* **Cancer_recapitulation_microarray**: Analysis showing results across studies consistent with results from **GTEx_TCGA_Comparisons**. Includes code for table S7.
   * **breast_large_n_recapitulate.rmd**: Analysis of GSE15852, called breast-1.
   * **breast_large_n_recapitulate_2.rmd**: Analysis of GSE109169, called breast-2.
   * **lung_recapitulate.Rmd**: Analysis of GSE31210, called lung-1.
   * **lung_recapitulate_2.Rmd**: Analysis of GSE19188, called lung-2.
   * **esophagus_recapitulate.rmd**: Analysis of GSE161533, called esophagus.
   * **colon_recapitulate.rmd**: Analysis of GSE44076, called colon.
* **GTEx_age_sex_analysis**: Code for all age and sex related analyses. Code to generate Fig 2A, 2C, Fig 3A, 3B, 3C, 3D, 3F, Fig S3, and Fig S4.
   * **age_sex_analysis_code_v8.rmd**: Analysis of limma results and ImageJ colon transverse results. Used to generate Fig 2A, 2C, 3A, 3B, and 3C.
   * **adipose_sex_analysis.rmd**: Analyzing the results that came from ImageJ and adipose tissue expression. Used to generate Fig 3D and 3F.
   * **top_bottom_transverse.rmd**: The code used to select individuals to analyze in the transverese colon analysis using ImageJ based on their normalized *ADIPOQ* expression.
* **Setup Code**: Extra code not used in the original analysis, but it allows the creation of various files used i nthe analysis itself.
   * **set_up_code.rmd**: An RMD files generating several files used in tha analysis:
      * **Generate gtex counts RDA**: Makes an rda file (gtex-gene-counts-v8.rda) from GTEx count data which includes:
         * dat:  Raw count gtex data with genes as rows and samples as columns
         * gtab: A two column dataframe including ENSG gene names and symbols
         * stab: A dataframe including sample data for all samples in the count data frame.
      * **Generate rda files used for the colon individual selector**: This section prepares vst normalized versions of colon transverse, colon sigmoid, adipose subcutaneous, and adipose visceral (omentum) used for selecting individuals in the file GTEx_age_sex_analysis/top_bottom_transverse.rmd.
