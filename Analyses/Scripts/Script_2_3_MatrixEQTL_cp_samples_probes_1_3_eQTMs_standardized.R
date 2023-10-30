########################################################################
## Title: This Script runs MatrixeQTL - eQTMs - for common samples and probes in CVS-Placenta & standardized variables
## Date: 01/23
## Author: Linda Dieckmann
########################################################################

## Section: load package(s)
########################################################################
library(MatrixEQTL)


########################################################################
## eQTM mapping 
########################################################################
base_dir <- getwd()

## Section: Parameter Settings
########################################################################
# Associations significant at pvOutputThreshold (pvOutputThreshold.cis) levels are saved to output_file_name (output_file_name.cis), 
# with corresponding estimates of effect size (slope coefficient), test statistics, p-values, and q-values (false discovery rate).
cis_threshold <- 5e-2  # cis-eQTMS cutoff: 0.05 = 5e-2; 0.01 = 1e-2; 0.00001 = 1e-5
cis_dist <- 1.5e5 #1.5e5 # Distance for local (cis) gene-SNP pairs: cis window of 1Mb = 1e6, 100kb = 1e5, 150kb = 1.5e5

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
expression_file_name <-  "02_Data/prepared/scaled_inv_norm_placenta_tmm_filtered_ordered_eqtm_common_cp_and_probes.txt"
methylation_file_name <- "02_Data/prepared/scaled_methylation_M_placenta_filtered_ordered_eqtm_common_cp_and_probes.txt"
covariates_file_name <- "02_Data/prepared/scaled_cov_placenta_ordered_eqtm_common_cp.txt"
  
gene_location_file_name <- "02_Data/prepared/gene_location_placenta_rna_common_probes.txt";
cpg_location_file_name <- "02_Data/prepared/cpg_location_placenta_meth_eqtm_common_probes.txt";

methylation_file_path <- file.path(base_dir, methylation_file_name)
expression_file_path <- file.path(base_dir, expression_file_name)
covariates_file_path <- file.path(base_dir, covariates_file_name)
gene_location_file_path <- file.path(base_dir, gene_location_file_name)
cpg_location_file_path <- file.path(base_dir, cpg_location_file_name)

output_file_name_cis = "02_Data/MatrixEQTL_Output/scaled_eqtm_placenta_cis_cp_samples_probes.txt"

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
scaled_me_eqtm_placenta_cp_samples_probes = Matrix_eQTL_main(
  snps = meth, # SlicedData object with methylation information. 
  gene = gene, # SlicedData object with gene expression information. 
  cvrt = cvrt, # SlicedData object with additional covariates.
  snpspos = cpgpos, # must have 3 columns - SNP name, chromosome, and position
  genepos = genepos, # data.frame with information about transcript locations, must have 4 columns - the name, chromosome, and positions of the left and right ends
  output_file_name.cis = output_file_name_cis, # local associations are saved to this file
  #output_file_name = output_file_name_tra, # significant associations are saved to this file (all significant associations if pvOutputThreshold=0 or only distant if pvOutputThreshold>0)
  pvOutputThreshold.cis = cis_threshold, # cis threshold for reporting
  pvOutputThreshold = 0, # trans threshold or 0 for only cis
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = FALSE,
  cisDist = cis_dist, # SNP-gene pairs within this distance are considered local. The distance is measured from the nearest end of the gene. SNPs within a gene are always considered local.
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE, # TRUE to record the minimum p-value for each SNP and each gene in the returned object. The minimum p-values are recorded even if if they are above the corresponding thresholds of pvOutputThreshold and pvOutputThreshold.cis. The analysis runs faster when the parameter is set to FALSE.
  noFDRsaveMemory = FALSE) # TRUE to save significant gene-SNP pairs directly to the output files, reduce memory footprint and skip FDR calculation. The eQTMs are not recorded in the returned object if noFDRsaveMemory = TRUE.

head(scaled_me_eqtm_placenta_cp_samples_probes$cis$eqtls)

save(scaled_me_eqtm_placenta_cp_samples_probes, file = "02_Data/MatrixEQTL_Output/scaled_me_eqtm_placenta_cp_samples_probes.Rdata")

## Section: delete placenta objects not needed anymore
########################################################################
rm(meth, gene, cvrt, cpgpos, genepos, scaled_me_eqtm_placenta_cp_samples_probes)

########################################################################
## Section: CVS
########################################################################

## Section: define file names and paths
########################################################################
methylation_file_name <-  "02_Data/prepared/scaled_methylation_M_cvs_filtered_ordered_eqtm_common_cp_and_probes.txt"
expression_file_name <- "02_Data/prepared/scaled_inv_norm_cvs_tmm_filtered_ordered_eqtm_common_cp_and_probes.txt"
covariates_file_name <- "02_Data/prepared/scaled_cov_cvs_ordered_eqtm_common_cp.txt"

cpg_location_file_name <- "02_Data/prepared/cpg_location_cvs_meth_eqtm_common_probes.txt";
gene_location_file_name <- "02_Data/prepared/gene_location_cvs_rna_common_probes.txt";

methylation_file_path <- file.path(base_dir, methylation_file_name)
expression_file_path <- file.path(base_dir, expression_file_name)
covariates_file_path <- file.path(base_dir, covariates_file_name)
cpg_location_file_path <- file.path(base_dir, cpg_location_file_name)
gene_location_file_path <- file.path(base_dir, gene_location_file_name)

output_file_name_cis = "02_Data/MatrixEQTL_Output/scaled_eqtm_cvs_cis_cp_samples_probes.txt"

## Section: load data
########################################################################
## Load meth data
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
scaled_me_eqtm_cvs_cp_samples_probes = Matrix_eQTL_main(
  snps = meth, # SlicedData object with genotype information.
  gene = gene, # SlicedData object with gene expression information.
  cvrt = cvrt, # SlicedData object with additional covariates.
  snpspos = cpgpos, # must have 3 columns - SNP name, chromosome, and position
  genepos = genepos, # data.frame with information about transcript locations, must have 4 columns - the name, chromosome, and positions of the left and right ends
  output_file_name.cis = output_file_name_cis, # local associations are saved to this file
  pvOutputThreshold.cis = cis_threshold,
  pvOutputThreshold = 0,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = FALSE,
  cisDist = cis_dist, # SNP-gene pairs within this distance are considered local. The distance is measured from the nearest end of the gene. SNPs within a gene are always considered local.
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE, # TRUE to record the minimum p-value for each SNP and each gene in the returned object. The minimum p-values are recorded even if if they are above the corresponding thresholds of pvOutputThreshold and pvOutputThreshold.cis. The analysis runs faster when the parameter is set to FALSE.
  noFDRsaveMemory = FALSE) # TRUE to save significant gene-SNP pairs directly to the output files, reduce memory footprint and skip FDR calculation. The eQTMs are not recorded in the returned object if noFDRsaveMemory = TRUE.

head(scaled_me_eqtm_cvs_cp_samples_probes$cis$eqtls)

save(scaled_me_eqtm_cvs_cp_samples_probes, file = "02_Data/MatrixEQTL_Output/scaled_me_eqtm_cvs_cp_samples_probes.Rdata")


