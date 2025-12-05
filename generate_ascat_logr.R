library(ASCAT)
library(argparse)

parser <- ArgumentParser(description = "Run ascat.prepareHTS on selected samples")
parser$add_argument("--sample_name_hash", help = "Name of the tumour sample, with sequencing hash")
parser$add_argument("--gl_hash",          help = "Name of the associated GL sample, with sequencing hash")
parser$add_argument("--sex",              help = "Sex of the patient")
parser$add_argument("--bam_path",         help = "Path to all bams")
parser$add_argument("--results_dir",      help = "Directory to store results")
parser$add_argument("--alleles_prefix",   help = "Prefix for the alleles files")
parser$add_argument("--loci_prefix",      help = "Prefix for the loci files")
parser$add_argument("--alignment_ext",    default = 'bam', help = "Alignment file extension, bam or cram")

args <-parser$parse_args()

# Generate inputs and outputs based on the sample name
tumour_bam <- paste0(args$bam_path, '/', args$sample_name_hash, '.', args$alignment_ext)
gl_bam <- paste0(args$bam_path, '/', args$gl_hash, '.', args$alignment_ext)

# Make sure BAMs exist
assertthat::assert_that(all(c(file.exists(tumour_bam), file.exists(gl_bam))), msg = 'Tumour and/or germline BAM/CRAMs are missing, check the BAM/CRAM path provided')
if (args$sex == 'Male') {
    clinical_sex <- 'XY'
} else if (args$sex == 'Female') {
    clinical_sex <- 'XX'
} else {
    stop('Argument --sex needs to be Male or Female')
}

# Output files
logr_tumour_file <- paste0(args$results_dir, '/', args$sample_name_hash, '_logr.txt')
logr_gl_file <- paste0(args$results_dir, '/', args$gl_hash, '_logr.txt')

baf_tumour_file <- paste0(args$results_dir, '/', args$sample_name_hash, '_baf.txt')
baf_gl_file <- paste0(args$results_dir, '/', args$gl_hash, '_baf.txt')

# Script to test the logR and BAF extraction from ASCAT
ascat.prepareHTS(
    tumourseqfile = tumour_bam,
    normalseqfile = gl_bam,
    tumourname = args$sample_name_hash,
    normalname = args$gl_hash,
    allelecounter_exe = "/camp/apps/eb/software/alleleCount/4.0.0-foss-2016b/bin/alleleCounter",
    alleles.prefix = args$alleles_prefix,
    loci.prefix = args$loci_prefix,
    gender = clinical_sex,
    genomeVersion = "hg38",
    nthreads = 2,
    tumourLogR_file = logr_tumour_file,
    tumourBAF_file = baf_tumour_file,
    normalLogR_file = logr_gl_file,
    normalBAF_file = baf_gl_file,
    minCounts = 20,
    BED_file = '/nemo/project/proj-tracerx-lung/tracerx/_PIPELINE/TRACERx-assets/v3/capture_targets/SureSelectV5/SureSelectv5_TRACERx_edition_hg38.padded.reduced.bed.gz',
    probloci_file = '/nemo/lab/swantonc/working/ruizc/2025-06-27-compare_hg38_qc/results/2025-10-30-test_ascat_logr/input/probloci.txt.gz',
    chrom_names = c(1:22, "X"),
    min_base_qual = 20,
    min_map_qual = 20,
    additional_allelecounter_flags = "--ref /nemo/project/proj-tracerx-lung/tracerx/_PIPELINE/TRACERx-assets/v3/reference/hg38/hg38_u2af1/GRCh38.d1.vd1.chr21fix.fa",
    skip_allele_counting_tumour = FALSE,
    skip_allele_counting_normal = FALSE,
    loci_binsize = 1,
    seed = 1234
)
