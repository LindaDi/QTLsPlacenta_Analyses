#!/bin/bash
#
#SBATCH --job-name="plink_clump_cvs_eqtl"
#SBATCH --ntasks=1
#SBATCH --array=1-12
#SBATCH --output=/binder/Linda/Processing/P2_Omics_PlacentaTissue/03_Code_Scripts/log_files/clumps_eqtl_cvs/array_%A-%a.log

# get all input data
FILES=(/binder/Linda/Processing/P2_Omics_PlacentaTissue/02_Data/clumping/cvs-eqtl_filesforclump/clump*.txt)    
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
# create basename of input file:
file_base=${basename ${FILE} .txt)

plink --bfile /binder/Linda/Data/ITU/genotypes/QCed_imputed_fullqced_maf_filter/itu_all_qc2_mafiltered --clump $FILE --clump-r2 0.2 --clump-kb 1000 --clump-p1 0.05 --clump-p2 1 --out /binder/Linda/Processing/P2_Omics_PlacentaTissue/02_Data/clumping/cvs-eqtl_filesfromclump/${file_base}



