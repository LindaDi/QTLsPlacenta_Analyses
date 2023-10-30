########################################################################
## Title: This Script runs BootstrapQTL - meQTLs - for common samples and probes in CVS-Placenta & standardized variables
## Date: 01/23
## Author: Linda Dieckmann
########################################################################

## Section: load package(s)
########################################################################
library(MatrixEQTL)
library(BootstrapQTL)

########################################################################
## meqtl mapping 
########################################################################
base_dir <- getwd()

## Section: Parameter Settings
########################################################################
# Associations significant at pvOutputThreshold (pvOutputThreshold.cis) levels are saved to output_file_name (output_file_name.cis), 
# with corresponding estimates of effect size (slope coefficient), test statistics, p-values, and q-values (false discovery rate).
cis_threshold <- 5e-2  # cis-meqtlS cutoff: 0.05 = 5e-2; 0.01 = 1e-2; 0.00001 = 1e-5
cis_dist <- 1.5e5 # Distance for local (cis) gene-SNP pairs: cis window of 1Mb = 1e6, 100kb = 1e5, 150kb = 1.5e5

# definition of parameters: 
  # Set pvOutputThreshold > 0 and pvOutputThreshold.cis = 0 (or use Matrix_eqtl_engine) to perform meqtl analysis without using gene/SNP locations. 
  # Associations significant at the pvOutputThreshold level are be recorded in output_file_name and in the returned object.
  # Set pvOutputThreshold = 0 and pvOutputThreshold.cis > 0 to perform meqtl analysis for local gene-SNP pairs only. 
  # Local associations significant at pvOutputThreshold.cis level will be recorded in output_file_name.cis and in the returned object.
  # Set pvOutputThreshold > 0 and pvOutputThreshold.cis > 0 to perform meqtl analysis with separate p-value thresholds for local and distant meqtls. 
  # Distant and local associations significant at corresponding thresholds are recorded in output_file_name and output_file_name.cis respectively and in the returned object. In this case the false discovery rate is calculated separately for these two sets of meqtls.

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR;
# Set useModel = modelLINEAR to model the effect of the genotype as additive linear and test for its significance using t-statistic.
# Set useModel = modelANOVA to treat genotype as a categorical variables and use ANOVA model and test for its significance using F-test. The default number of ANOVA categories is 3. Set otherwise like this: options(MatrixEQTL.ANOVA.categories=4).
# Set useModel = modelLINEAR_CROSS to add a new term to the modelequal to the product of genotype and the last covariate; the significance of this term is then tested using t-statistic.

# Error covariance matrix
errorCovariance = numeric()
# Set to numeric() for identity

#######################################################################
# Section: Placenta
#######################################################################

## Section: define file names and paths
########################################################################
SNP_file_name <-  "../../02_Data/prepared/scaled_geno_t_fullqced_IDc_placenta_ordered_meqtl_common_cp.txt"
methylation_file_name <- "../../02_Data/prepared/scaled_methylation_M_placenta_filtered_ordered_meqtl_common_cp_and_probes.txt"
covariates_file_name <- "../../02_Data/prepared/scaled_cov_meth_placenta_ordered_meqtl_common_cp.txt"

snps_location_file_name <- "../../02_Data/prepared/snplocation.txt";
cpg_location_file_name <- "../../02_Data/prepared/cpg_location_placenta_meth_meqtl_common_probes.txt";

SNP_file_path <- file.path(base_dir, SNP_file_name)
methylation_file_path <- file.path(base_dir, methylation_file_name)
covariates_file_path <- file.path(base_dir, covariates_file_name)
snps_location_file_path <- file.path(base_dir, snps_location_file_name)
cpg_location_file_path <- file.path(base_dir, cpg_location_file_name)

output_file_name_cis = "scaled_meqtl_placenta_cis_cp_samples_probes.txt"

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

## Load gene methylation data
gene = SlicedData$new()
gene$fileDelimiter = "\t" # the TAB character
gene$fileOmitCharacters = "NA" # denote missing values;
gene$fileSkipRows = 1 # one row of column labels
gene$fileSkipColumns = 1 # one column of row labels
gene$fileSliceSize = 2000 # read file in slices of 2,000 rows
gene$LoadFile(methylation_file_name)

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
cpgpos = read.table(cpg_location_file_name, header = TRUE, stringsAsFactors = FALSE)

## Section: run analysis
########################################################################
# Run the BootstrapQTL analysis
scaled_bme_meqtl_placenta_cp_samples_probes <- BootstrapQTL(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  snpspos = snpspos,
  genepos = cpgpos,
  n_bootstraps=200, n_cores=10, # default 200 bootstraps
  eGene_detection_file_name = output_file_name_cis, # File to save local cis associations to in the eGene detection analysis. Corresponds to the output_file_name.cis argument in Matrix_eQTL_main.
  bootstrap_file_directory = "../../02_Data/MatrixEQTL_Output/bootstrap_analyses/meQTL_placenta/",
  useModel = useModel,
  errorCovariance = errorCovariance,
  cisDist = cis_dist,
  local_correction = "bonferroni", # is standard, another alternative would be eigenMT -> then 'eigenMT_tests_per_gene' must be a data.frame containing the number of effective independent tests per gene
  global_correction = "BH",
  correction_type = "shrinkage" # default
)

save(scaled_bme_meqtl_placenta_cp_samples_probes, file = "../../02_Data/MatrixEQTL_Output/scaled_bme_meqtl_placenta_cp_samples_probes.Rdata")

# Section: delete placenta objects not needed anymore
#######################################################################
rm(snps, gene, cvrt, snpspos, cpgpos, scaled_bme_meqtl_placenta_cp_samples_probes)

########################################################################
## Section: CVS
########################################################################

## Section: define file names and paths
########################################################################
SNP_file_name <- "../../02_Data/prepared/scaled_geno_t_fullqced_IDc_cvs_ordered_meqtl_common_cp.txt"
methylation_file_name <- "../../02_Data/prepared/scaled_methylation_M_cvs_filtered_ordered_meqtl_common_cp_and_probes.txt"
covariates_file_name <- "../../02_Data/prepared/scaled_cov_meth_cvs_ordered_meqtl_common_cp.txt"

snps_location_file_name <- "../../02_Data/prepared/snplocation.txt";
cpg_location_file_name <- "../../02_Data/prepared/cpg_location_cvs_meth_meqtl_common_probes.txt";

SNP_file_path <- file.path(base_dir, SNP_file_name)
methylation_file_path <- file.path(base_dir, methylation_file_name)
covariates_file_path <- file.path(base_dir, covariates_file_name)
snps_location_file_path <- file.path(base_dir, snps_location_file_name)
cpg_location_file_path <- file.path(base_dir, cpg_location_file_name)

output_file_name_cis = "scaled_meqtl_cvs_cis_cp_samples_probes.txt"

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

## Load gene methylation data
gene = SlicedData$new()
gene$fileDelimiter = "\t" # the TAB character
gene$fileOmitCharacters = "NA" # denote missing values;
gene$fileSkipRows = 1 # one row of column labels
gene$fileSkipColumns = 1 # one column of row labels
gene$fileSliceSize = 2000 # read file in slices of 2,000 rows
gene$LoadFile(methylation_file_name)

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
cpgpos = read.table(cpg_location_file_name, header = TRUE, stringsAsFactors = FALSE)

## Section: run analysis
########################################################################
# Run the BootstrapQTL analysis
scaled_bme_meqtl_cvs_cp_samples_probes <- BootstrapQTL(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  snpspos = snpspos,
  genepos = cpgpos,
  n_bootstraps=200, n_cores=10, # default 200 bootstraps
  eGene_detection_file_name = output_file_name_cis, # File to save local cis associations to in the eGene detection analysis. Corresponds to the output_file_name.cis argument in Matrix_eQTL_main.
  bootstrap_file_directory = "../../02_Data/MatrixEQTL_Output/bootstrap_analyses/meQTL_cvs/",
  useModel = useModel,
  errorCovariance = errorCovariance,
  cisDist = cis_dist,
  local_correction = "bonferroni", # is standard, another alternative would be eigenMT -> then 'eigenMT_tests_per_gene' must be a data.frame containing the number of effective independent tests per gene
  global_correction = "BH",
  correction_type = "shrinkage" # default
)

save(scaled_bme_meqtl_cvs_cp_samples_probes, file = "../../02_Data/MatrixEQTL_Output/scaled_bme_meqtl_cvs_cp_samples_probes.Rdata")

########################################################################




