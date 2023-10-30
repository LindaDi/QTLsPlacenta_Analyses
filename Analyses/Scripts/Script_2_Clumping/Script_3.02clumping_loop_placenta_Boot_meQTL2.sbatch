#!/bin/bash
#SBATCH --output=/binder/Linda/Processing/P2_Omics_PlacentaTissue/03_Code_Scripts/log_files/clumps_meqtl_placenta/1.2/loop-submission_%A-%a.log
#SBATCH --job-name=loop_scripts

rslt_dir=/home/ldieckmann/Processing/P2_Omics_PlacentaTissue/placenta_meqtl_clump_output_2

for i in 0 10000 20000
do
   sbatch ./Script_3.02clumping_placenta_Boot_meQTL_get_snps_per_cpg2.sh $i
done



