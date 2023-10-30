########################################################################
## Title: This Script runs BootstrapQTL - for Cluster - eQTMs (expression, methylation)
## Date: 10/12/21
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

########################################################################
## eQTM mapping 
########################################################################
base_dir <- getwd()

## Section: Parameter Settings
########################################################################
# Associations significant at pvOutputThreshold (pvOutputThreshold.cis) levels are saved to output_file_name (output_file_name.cis), 
# with corresponding estimates of effect size (slope coefficient), test statistics, p-values, and q-values (false discovery rate).
cis_threshold <- 5e-2  # cis-eQTMS cutoff: 0.05 = 5e-2; 0.01 = 1e-2; 0.00001 = 1e-5
cis_dist <- 1.5e5 # Distance for local (cis) gene-SNP pairs: cis window of 1Mb = 1e6, 100kb = 1e5, 150kb = 1.5e5


# definition of parameters: 
  # Set pvOutputThreshold > 0 and pvOutputThreshold.cis = 0 (or use Matrix_eQTL_engine) to perform eQTM analysis without using gene/SNP locations. 
  # Associations significant at the pvOutputThreshold level are be recorded in output_file_name and in the returned object.
  # Set pvOutputThreshold = 0 and pvOutputThreshold.cis > 0 to perform eQTM analysis for local gene-SNP pairs only. 
  # Local associations significant at pvOutputThreshold.cis level will be recorded in output_file_name.cis and in the returned object.
  # Set pvOutputThreshold > 0 and pvOutputThreshold.cis > 0 to perform eQTM analysis with separate p-value thresholds for local and distant eQTMs. 
  # Distant and local associations significant at corresponding thresholds are recorded in output_file_name and output_file_name.cis respectively and in the returned object. In this case the false discovery rate is calculated separately for these two sets of eQTMs.

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR;
# Set useModel = modelLINEAR to model the effect of the genotype as additive linear and test for its significance using t-statistic.
# Set useModel = modelANOVA to treat genotype as a categorical variables and use ANOVA model and test for its significance using F-test. The default number of ANOVA categories is 3. Set otherwise like this: options(MatrixeQTM.ANOVA.categories=4).
# Set useModel = modelLINEAR_CROSS to add a new term to the modelequal to the product of genotype and the last covariate; the significance of this term is then tested using t-statistic.

# Error covariance matrix
errorCovariance = numeric()
# Set to numeric() for identity

########################################################################
## Section: Placenta
########################################################################

## Section: define file names and paths
########################################################################
expression_file_name <-  "../02_Data/prepared/inv_norm_placenta_tmm_filtered_ordered_eqtm.txt";
methylation_file_name <- "../02_Data/prepared/methylation_M_placenta_filtered_ordered_eqtm.txt";
covariates_file_name <- "../02_Data/prepared/cov_placenta_ordered_eqtm.txt";
gene_location_file_name <- "../02_Data/prepared/gene_location_placenta_rna.txt";
cpg_location_file_name <- "../02_Data/prepared/cpg_location_placenta_meth_eqtm.txt";

methylation_file_path <- file.path(base_dir, methylation_file_name)
expression_file_path <- file.path(base_dir, expression_file_name)
covariates_file_path <- file.path(base_dir, covariates_file_name)
gene_location_file_path <- file.path(base_dir, gene_location_file_name)
cpg_location_file_path <- file.path(base_dir, cpg_location_file_name)

output_file_name_cis = "../02_Data/MatrixEQTL_Output/eqtm_placenta_cis.txt"

## Section: load data
########################################################################
## Load methylation data
meth <- SlicedData$new()
meth$fileDelimiter = "\t" # the TAB character
meth$fileOmitCharacters = "NA" # denote missing values
meth$fileSkipRows = 1 # one row of column labels
meth$fileSkipColumns = 1 # one column of row labels
meth$fileSliceSize = 2000 # read file in slices of 2,000 rows
meth$LoadFile(methylation_file_name)

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

cpgpos = read.table(cpg_location_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)

## Section: run analysis
########################################################################
# Run the BootstrapQTL analysis
bme_eqtm_placenta <- BootstrapQTL(
  snps = meth, 
  gene = gene, 
  cvrt = cvrt,
  snpspos = cpgpos, 
  genepos = genepos,
  n_bootstraps=200, n_cores=10, # default 200 bootstraps
  eGene_detection_file_name = output_file_name_cis, # File to save local cis associations to in the eGene detection analysis. Corresponds to the output_file_name.cis argument in Matrix_eQTL_main. 
  bootstrap_file_directory = "../02_Data/MatrixEQTL_Output/bootstrap_analyses/eQTM_placenta/",
  useModel = useModel,
  errorCovariance = errorCovariance,
  cisDist = cis_dist,
  local_correction = "bonferroni", # is standard, another alternative would be eigenMT -> then 'eigenMT_tests_per_gene' must be a data.frame containing the number of effective independent tests per gene
  global_correction = "BH", # default
  correction_type = "shrinkage" # default
)

save(bme_eqtm_placenta, file = "../02_Data/MatrixEQTL_Output/bme_eqtm_placenta.Rdata")


## Section: delete placenta objects not needed anymore
########################################################################
rm(meth, gene, cvrt, cpgpos, genepos, me_eqtm_placenta)

########################################################################
## Section: CVS
########################################################################

## Section: define file names and paths
########################################################################
methylation_file_name <-  "../02_Data/prepared/methylation_M_cvs_filtered_ordered_eqtm.txt";
expression_file_name <- "../02_Data/prepared/inv_norm_cvs_tmm_filtered_ordered_eqtm.txt";
covariates_file_name <- "../02_Data/prepared/cov_cvs_ordered_eqtm.txt";
cpg_location_file_name <- "../02_Data/prepared/cpg_location_cvs_meth_eqtm.txt";
gene_location_file_name <- "../02_Data/prepared/gene_location_cvs_rna.txt";

methylation_file_path <- file.path(base_dir, methylation_file_name)
expression_file_path <- file.path(base_dir, expression_file_name)
covariates_file_path <- file.path(base_dir, covariates_file_name)
cpg_location_file_path <- file.path(base_dir, cpg_location_file_name)
gene_location_file_path <- file.path(base_dir, gene_location_file_name)

output_file_name_cis = "eqtm_cvs_cis.txt"

## Section: load data
########################################################################
## Load genotype data
meth <- SlicedData$new()
meth$fileDelimiter = "\t" # the TAB character
meth$fileOmitCharacters = "NA" # denote missing values
meth$fileSkipRows = 1 # one row of column labels
meth$fileSkipColumns = 1 # one column of row labels
meth$fileSliceSize = 2000 # read file in slices of 2,000 rows
meth$LoadFile(methylation_file_name)

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

cpgpos = read.table(cpg_location_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)

## Section: run analysis
########################################################################
# Run the BootstrapQTL analysis
bme_eqtm_cvs <- BootstrapQTL(
  snps = meth, 
  gene = gene, 
  cvrt = cvrt,
  snpspos = cpgpos, 
  genepos = genepos,
  n_bootstraps=200, n_cores=10, # default 200 bootstraps
  eGene_detection_file_name = output_file_name_cis, # File to save local cis associations to in the eGene detection analysis. Corresponds to the output_file_name.cis argument in Matrix_eQTL_main. 
  bootstrap_file_directory = "../02_Data/MatrixEQTL_Output/bootstrap_analyses/eQTM_cvs/",
  useModel = useModel,
  errorCovariance = errorCovariance,
  cisDist = cis_dist,
  local_correction = "bonferroni", # is standard, another alternative would be eigenMT -> then 'eigenMT_tests_per_gene' must be a data.frame containing the number of effective independent tests per gene
  global_correction = "BH", # default
  correction_type = "shrinkage" # default
)

save(bme_eqtm_cvs, file = "../02_Data/MatrixEQTL_Output/bme_eqtm_cvs.Rdata")


