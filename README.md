# Matrisome_GTEx_analysis
The Github Repository for the paper "Sex, age, tissue, and disease patterns of matrisome expression in GTEx transcriptome data" by Tim O. Nieuwenhuis, Avi Z. Rosenberg, Matthew N. McCall, and Marc K. Halushka. 

## Purpose
This Github Repo exists for the purpose of giving individuals the ability to reproduce the findings and figures found in the previously named paper. All code was written by Tim Nieuwenhuis and as a result any issues with the Github should be directed towards him via issues or email correspodence. 

## Organization
* **global_in**: Holds files used in multiple analyses
    *  matrisome_hs_masterlist_EDIT_r.csv: A csv version of the list of Matrisome Genes made by the Hynes lab. Original data can be retrieved [here](http://matrisomeproject.mit.edu/static/media/uploads/Logos/.thumbnails/excel-2010-icon.jpg/excel-2010-icon-40x30.jpg).
    *  **tim_scripts.R**: A file holding small quality of life scripts that are used when working with GTEx data.
    *  **gtex_phenotypes_v8.csv** *not available*: This is a file that includes specific phenotype data related to the individuals in GTEx such as their true age. This data is not included in the git hub as it is only accessible through DBGap with the correct permissions. The main usage of this tool is for non-binned ages of the individuals in GTEx, however one can use the binned ages available [here](https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.xlsx) as a substitute. 
* **GTEx_protein_gene_correlations**: Looking at GTEx Specific protein gene correlations and correlations from Wang et al. Includes code for Figure 1A, Figure 1B, Figure S1, and Figure S2.
  * **protein_gene_correlation_analysis**: This file includes the analysis of correlations from [Jiang et al.](https://www.sciencedirect.com/science/article/pii/S0092867420310783). The associated data files from the paper required to run the code are to large to upload to Github, and so their https are found in the comments. This file includes code for Figure S1.
  * ADD SPECIFIC FILE HERE
* **GTEx_intra_tissue_comparisons**: 
* **GTEx_median_tissue_comparisons**:
* **GTEx_TCGA_comparison**:
* **GTEx_IPF_analysis**:
