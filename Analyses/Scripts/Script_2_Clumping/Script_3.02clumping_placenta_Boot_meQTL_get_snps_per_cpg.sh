#!/bin/bash
#SBATCH --job-name=c_meqtl_p
#SBATCH --output=/binder/Linda/Processing/P2_Omics_PlacentaTissue/03_Code_Scripts/log_files/clumps_meqtl_placenta/2/ind_submission_%A-%a.logs
#SBATCH --mem-per-cpu=200M     # Each task uses max 200M of memory
#SBATCH --array=1-10000%100    # Submit 10'000 tasks. Run max 100 concurrently
#SBATCH --part=pe
#SBATCH --exclude=pe5,pe9,pe8

snp_dir=/binder/Linda/Data/ITU/genotypes/QCed_imputed_fullqced_maf_filter
input_dir=/binder/Linda/Processing/P2_Omics_PlacentaTissue/02_Data/clumping/placenta-meqtl_filesforclump
rslt_dir=/home/ldieckmann/Processing/P2_Omics_PlacentaTissue/placenta_meqtl_clump_output

SCRATCH_DIRECTORY=${rslt_dir}/scratch/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}

i=$1
n=$(($SLURM_ARRAY_TASK_ID+${i}))
cpg=$(awk -v cpg_id=$n 'NR==cpg_id {print $1}' ${input_dir}/unique_eCpGs_placenta_for_clump.txt)

echo "SNP P" > snp_lst_for_ld_clump_${cpg}.txt
awk -v cg=$cpg '$1 ~ cg {print $2 " " $9}' ${input_dir}/bme_meqtl_placenta_for_clump.txt >> snp_lst_for_ld_clump_${cpg}.txt

plink --bfile $snp_dir/itu_all_qc2_mafiltered --clump snp_lst_for_ld_clump_${cpg}.txt --clump-r2 0.2 --clump-kb 150 --clump-p1 0.05 --clump-p2 1 --out clump_${cpg}

awk -v cg=$cpg 'NR != 1 {if($3) print cg " " $3}' clump_${cpg}.clumped >> ${SCRATCH_DIRECTORY}/me-qtl_cis_ind_cpg_snp_associations_${cpg}.txt

cp ${SCRATCH_DIRECTORY}/me-qtl_cis_ind_cpg_snp_associations_${cpg}.txt ${rslt_dir}

cd ${rslt_dir}
rm -rf ${SCRATCH_DIRECTORY}

exit 0

