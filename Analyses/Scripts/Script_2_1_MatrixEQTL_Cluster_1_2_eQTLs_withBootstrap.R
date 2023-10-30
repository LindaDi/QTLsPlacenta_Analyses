########################################################################
## Title: This Script runs BootstrapQTL - for Cluster - eQTLs
## Date: 14/01/22
## Author: Linda Dieckmann
########################################################################

## Section: load package(s)
########################################################################
library(MatrixEQTL)
library(BootstrapQTL)
library(data.table)
library(parallel)
library(doParallel)
library(doMC)
# documentation: https://cran.r-project.org/web/packages/BootstrapQTL/BootstrapQTL.pdf

########################################################################
## eQTL mapping 
########################################################################
base_dir <- getwd()

## Section: Parameter Settings
########################################################################
# Associations significant at pvOutputThreshold (pvOutputThreshold.cis) levels are saved to output_file_name (output_file_name.cis), 
# with corresponding estimates of effect size (slope coefficient), test statistics, p-values, and q-values (false discovery rate).
cis_threshold <- 5e-2  # cis-eQTLS cutoff: 0.05 = 5e-2; 0.01 = 1e-2; 0.00001 = 1e-5
# trans_threshold <- 1e-5 # trans-eQTLs cutoff: 0.05 = 5e-2; 0.01 = 1e-2 (0 means no trans-eQTLs); 0.00005 = 5e-5; 0.00001 = 1e-5
cis_dist <- 1e6 # Distance for local (cis) gene-SNP pairs: cis window of 1Mb = 1e6, 100kb = 1e5


# definition of parameters: 
  # Set pvOutputThreshold > 0 and pvOutputThreshold.cis = 0 (or use Matrix_eQTL_engine) to perform eQTL analysis without using gene/SNP locations. 
  # Associations significant at the pvOutputThreshold level are be recorded in output_file_name and in the returned object.
  # Set pvOutputThreshold = 0 and pvOutputThreshold.cis > 0 to perform eQTL analysis for local gene-SNP pairs only. 
  # Local associations significant at pvOutputThreshold.cis level will be recorded in output_file_name.cis and in the returned object.
  # Set pvOutputThreshold > 0 and pvOutputThreshold.cis > 0 to perform eQTL analysis with separate p-value thresholds for local and distant eQTLs. 
  # Distant and local associations significant at corresponding thresholds are recorded in output_file_name and output_file_name.cis respectively and in the returned object. In this case the false discovery rate is calculated separately for these two sets of eQTLs.

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR;
# Set useModel = modelLINEAR to model the effect of the genotype as additive linear and test for its significance using t-statistic.
# Set useModel = modelANOVA to treat genotype as a categorical variables and use ANOVA model and test for its significance using F-test. The default number of ANOVA categories is 3. Set otherwise like this: options(MatrixEQTL.ANOVA.categories=4).
# Set useModel = modelLINEAR_CROSS to add a new term to the modelequal to the product of genotype and the last covariate; the significance of this term is then tested using t-statistic.

# Error covariance matrix
errorCovariance = numeric()
# Set to numeric() for identity

########################################################################
## Section: Placenta
########################################################################

## Section: define file names and paths
########################################################################
SNP_file_name <-  "../02_Data/prepared/geno_t_fullqced_IDc_placenta_ordered_eqtl.txt";
expression_file_name <- "../02_Data/prepared/inv_norm_placenta_tmm_filtered_ordered_eqtl.txt";
covariates_file_name <- "../02_Data/prepared/cov_rna_placenta_ordered_eqtl.txt";
snps_location_file_name <- "../02_Data/prepared/snplocation.txt";
gene_location_file_name <- "../02_Data/prepared/gene_location_placenta_rna.txt";

SNP_file_path <- file.path(base_dir, SNP_file_name)
expression_file_path <- file.path(base_dir, expression_file_name)
covariates_file_path <- file.path(base_dir, covariates_file_name)
snps_location_file_path <- file.path(base_dir, snps_location_file_name)
gene_location_file_path <- file.path(base_dir, gene_location_file_name)

output_file_name_cis = "../02_Data/MatrixEQTL_Output/beqtl_placenta_cis.txt"

## Section: load data
########################################################################
## Load genotype data
snps <- SlicedData$new()
snps$fileDelimiter = "\t" # the TAB character
snps$fileOmitCharacters = "NA" # denote missing values
snps$fileSkipRows = 1 # one row of column labels
snps$fileSkipColumns = 1 # one column of row labels
snps$fileSliceSize = 2000 # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name)

## Load gene expression data
gene = SlicedData$new()
gene$fileDelimiter = "\t" # the TAB character
gene$fileOmitCharacters = "NA" # denote missing values;
gene$fileSkipRows = 1 # one row of column labels
gene$fileSkipColumns = 1 # one column of row labels
gene$fileSliceSize = 2000 # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name)

## Load covariates
cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"      # the TAB character
cvrt$fileOmitCharacters = "NA" # denote missing values;
cvrt$fileSkipRows = 1          # one row of column labels
cvrt$fileSkipColumns = 1      # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name)
}

snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)

## Section: run analysis
########################################################################
# Run the BootstrapQTL analysis
bme_eqtl_placenta <- BootstrapQTL(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  snpspos = snpspos, 
  genepos = genepos,
  n_bootstraps=200, n_cores=10, # default 200 bootstraps
  eGene_detection_file_name = output_file_name_cis, # File to save local cis associations to in the eGene detection analysis. Corresponds to the output_file_name.cis argument in Matrix_eQTL_main. 
  bootstrap_file_directory = "../02_Data/MatrixEQTL_Output/bootstrap_analyses/eQTL_placenta/",
  useModel = useModel,
  errorCovariance = errorCovariance,
  cisDist = cis_dist,
  local_correction = "bonferroni", # is standard, another alternative would be eigenMT -> then 'eigenMT_tests_per_gene' must be a data.frame containing the number of effective independent tests per gene
  global_correction = "BH", # default
  correction_type = "shrinkage" # default
)

save(bme_eqtl_placenta, file = "../02_Data/MatrixEQTL_Output/bme_eqtl_placenta.Rdata")

## Section: delete placenta objects not needed anymore
########################################################################
rm(snps, gene, cvrt, snpspos, genepos, bme_eqtl_placenta)


########################################################################
## Section: CVS
########################################################################

## Section: define file names and paths
########################################################################
SNP_file_name <-  "../02_Data/prepared/geno_t_fullqced_IDc_cvs_ordered_eqtl.txt";
expression_file_name <- "../02_Data/prepared/inv_norm_cvs_tmm_filtered_ordered_eqtl.txt";
covariates_file_name <- "../02_Data/prepared/cov_rna_cvs_ordered_eqtl.txt";
snps_location_file_name <- "../02_Data/prepared/snplocation.txt";
gene_location_file_name <- "../02_Data/prepared/gene_location_cvs_rna.txt";

SNP_file_path <- file.path(base_dir, SNP_file_name)
expression_file_path <- file.path(base_dir, expression_file_name)
covariates_file_path <- file.path(base_dir, covariates_file_name)
snps_location_file_path <- file.path(base_dir, snps_location_file_name)
gene_location_file_path <- file.path(base_dir, gene_location_file_name)

output_file_name_cis = "../02_Data/MatrixEQTL_Output/beqtl_cvs_cis.txt"

## Section: load data
########################################################################
## Load genotype data
snps <- SlicedData$new()
snps$fileDelimiter = "\t" # the TAB character
snps$fileOmitCharacters = "NA" # denote missing values
snps$fileSkipRows = 1 # one row of column labels
snps$fileSkipColumns = 1 # one column of row labels
snps$fileSliceSize = 2000 # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name)

## Load gene expression data
gene = SlicedData$new()
gene$fileDelimiter = "\t" # the TAB character
gene$fileOmitCharacters = "NA" # denote missing values;
gene$fileSkipRows = 1 # one row of column labels
gene$fileSkipColumns = 1 # one column of row labels
gene$fileSliceSize = 2000 # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name)

## Load covariates
cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"      # the TAB character
cvrt$fileOmitCharacters = "NA" # denote missing values;
cvrt$fileSkipRows = 1          # one row of column labels
cvrt$fileSkipColumns = 1      # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name)
}

snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)

## Section: run analysis
########################################################################
# Run the BootstrapQTL analysis
bme_eqtl_cvs <- BootstrapQTL(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  snpspos = snpspos, 
  genepos = genepos,
  n_bootstraps=200, n_cores=10, # default 200 bootstraps
  eGene_detection_file_name = output_file_name_cis, # File to save local cis associations to in the eGene detection analysis. Corresponds to the output_file_name.cis argument in Matrix_eQTL_main. 
  bootstrap_file_directory = "../02_Data/MatrixEQTL_Output/bootstrap_analyses/eQTL_cvs/",
  useModel = useModel,
  errorCovariance = errorCovariance,
  cisDist = cis_dist,
  local_correction = "bonferroni", # is standard, another alternative would be eigenMT -> then 'eigenMT_tests_per_gene' must be a data.frame containing the number of effective independent tests per gene
  global_correction = "BH", #default
  correction_type = "shrinkage" # default
)

save(bme_eqtl_cvs, file = "../02_Data/MatrixEQTL_Output/bme_eqtl_cvs.Rdata")



