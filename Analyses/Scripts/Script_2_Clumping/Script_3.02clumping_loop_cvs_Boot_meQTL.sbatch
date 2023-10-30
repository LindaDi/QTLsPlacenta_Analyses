#!/bin/bash
#SBATCH --output=/binder/Linda/Processing/P2_Omics_PlacentaTissue/03_Code_Scripts/log_files/clumps_meqtl_cvs/1/loop-submission_%A-%a.log
#SBATCH --job-name=loop_scripts

rslt_dir=/home/ldieckmann/Processing/P2_Omics_PlacentaTissue/cvs_meqtl_clump_output

for i in 0 10000 20000 30000 40000 50000
do
   sbatch ./Script_3.02clumping_cvs_Boot_meQTL_get_snps_per_cpg.sh $i
done



