#!/bin/bash
#SBATCH --job-name=ascat_hts             # Job name
#SBATCH --ntasks=1                              # Run one tasks
#SBATCH --mem=16G                               # Job Memory
#SBATCH --time=11:00:00                          # Time limit hrs:min:sec
#SBATCH --output=tmp/logs/ascat_hts_%A_%a.log     # Standard output and error log
#SBATCH --array=1-5                         #Number of tasks to run 

ml purge
module load R/4.0.0-foss-2020a
module load alleleCount/4.0.0-foss-2016b

# Get sample if interest, GL and sex
ROW=$(sed -n "${SLURM_ARRAY_TASK_ID}p" tmp/samples.txt)
SAMPLE_NAME=$(echo ${ROW} | cut -d ' ' -f 1)
PATIENT=$(echo ${SAMPLE_NAME} | cut -d '_' -f 1)
GL_NAME=$(echo ${ROW} | cut -d ' ' -f 2)
SEX=$(echo ${ROW} | cut -d ' ' -f 3)
VCF_PATH=/nemo/project/proj-tracerx-lung/tracerx/_RELEASE/release_tx842_hg38/${PATIENT}/germline/deepvariant/${GL_NAME}.deepvariant.vcf.gz
BAM_PATH=/nemo/project/proj-tracerx-lung/datasets/tracerx/alignments/wes/v3/bam
VARSCAN_PATH=/nemo/project/proj-tracerx-lung/tracerx/_RELEASE/release_tx842/${PATIENT}/ascat/${SAMPLE_NAME}/ascat_default/${SAMPLE_NAME}.ascat_default.ascat.bc.rds

# Generate results directory pper sample
RESULTS_OUT=$(echo ${PWD}/results/${SAMPLE_NAME})
mkdir -p ${RESULTS_OUT}/snp_input
mkdir -p ${RESULTS_OUT}/plot_comparisons

# Generate SNP files per sample
python generate_loci_allele_files.py --input_vcf ${VCF_PATH} \
                                     --loci_input input/G1000_allelesAll_hg38/ \
                                     --caller deepvariant \
                                     --out_dir ${RESULTS_OUT}/snp_input

cd ${RESULTS_OUT}
# Generate BAF and logR
Rscript ../../generate_ascat_logr.R --sample_name_hash ${SAMPLE_NAME} \
                              --gl_hash ${GL_NAME} \
                              --sex ${SEX} \
                              --bam_path ${BAM_PATH} \
                              --results_dir ${RESULTS_OUT} \
                              --alleles_prefix ${RESULTS_OUT}'/snp_input/alleles_chr' \
                              --loci_prefix ${RESULTS_OUT}'/snp_input/loci_chr'  \
                              --alignment_ext 'cram'

# Generate comparison plots against previous versions of BAF and logR
Rscript ../../test_logr_baf_diffs.R --varscan_obj_path ${VARSCAN_PATH} \
                              --sample_name_hash ${SAMPLE_NAME} \
                              --gl_name_hash ${GL_NAME} \
                              --ascat_baf ${RESULTS_OUT}/${SAMPLE_NAME}_baf.txt \
                              --ascat_logr ${RESULTS_OUT}/${SAMPLE_NAME}_logr.txt \
                              --results_dir ${RESULTS_OUT} \
                              --lift_hg38 TRUE \
                              --chain_file_path /flask/scratch/swantonc/ruizc/test_baf_logr/hg38ToHg19.over.chain

