library(data.table)
library(ggplot2)
library(GenomicRanges)

# Functions
filter_exon_capture <- function(to_filter_dt, exon_capture) {
    to_filter_gr <- GRanges(
        seqnames = to_filter_dt$chr,
        ranges = IRanges(start = to_filter_dt$pos, end = to_filter_dt$pos))
    # Overlap
    filtered_dt <- unique(as.data.table(mergeByOverlaps(to_filter_gr, exon_capture)))

    # Get columns of interest
    old_cols <- c('to_filter_gr.seqnames', 'to_filter_gr.start')
    filtered_dt <- copy(unique(filtered_dt[, ..old_cols]))
    setnames(filtered_dt, old_cols, c('chr', 'pos'))
    
    return(filtered_dt)
}


# Script to compare the 421 problematic loci file vs the "standard" Battenberg problematic loci file
probloci_421 <- fread('/nemo/project/proj-tracerx-lung/tracerx/_PIPELINE/TRACERx-assets/v1/ascat/probloci.tracerx.platypus.20210122.txt')
probloci_bat <-fread('input/probloci_270415.txt')
setnames(probloci_bat, c('Chr', 'Pos'), c('chr', 'pos'))
probloci_bat[, chr := gsub('X','23', chr)]

# Keep only positions in the exon capture kit
capture_kit <- fread('/nemo/project/proj-tracerx-lung/tracerx/_PIPELINE/TRACERx-assets/v1/capture_targets/SureSelectV5_TRACERx_Edition.padded.reduced.bed')
capture_kit_granges <- GRanges(
  seqnames = gsub('Y','24',gsub('X','23', gsub('chr', '', capture_kit$V1))),
  ranges = IRanges(start = capture_kit$V2, end = capture_kit$V3)
)

filt_probloci_421 <- filter_exon_capture(probloci_421, capture_kit_granges)
filt_probloci_421[, source := 'tx421']
filt_probloci_bat <- filter_exon_capture(probloci_bat, capture_kit_granges)
filt_probloci_bat[, source := 'battenberg']

# Find positions overlapping and missing in either
merged_probloci <- rbind(filt_probloci_421, filt_probloci_bat)

nb_probloci <- merged_probloci[, nb_calls := .N, by = c('chr', 'pos')]
nb_probloci[, called_by := ifelse(nb_calls == 2, 'both', source)]
nb_probloci <- unique(nb_probloci[, c('chr', 'pos', 'called_by')])
nb_probloci <- nb_probloci[, .(nb_calls = .N), by = c('chr', 'called_by')]

# Get total per caller removing chr Y (24), only present in 421 calls
nb_probloci_total <- nb_probloci[chr != 24, .(total_nb = sum(nb_calls)), by = called_by]

pdf('results/compare_probloci_snps.pdf')
ggplot(nb_probloci, aes(x = chr, y = nb_calls, fill = called_by)) +
    geom_bar(stat = 'identity', position = position_dodge()) +
    theme_bw() +
    ylab('# SNPs') +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size = 14, family = "ArialMT"),
          axis.text.y  = element_text(size = 12, family = "ArialMT"),
          axis.text.x  = element_text(size = 12, family = "ArialMT"),
          axis.title.x = element_text(size = 14, family = "ArialMT"),
          legend.title = element_text(size = 14, family = "ArialMT"),
          legend.text = element_text(size = 12, family = "ArialMT"),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 12, family = "ArialMT"))

ggplot(nb_probloci_total, aes(x = called_by, y = total_nb, fill = called_by)) +
    geom_bar(stat = 'identity', position = position_dodge()) +
    theme_bw() +
    ylab('# SNPs') + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size = 14, family = "ArialMT"),
          axis.text.y  = element_text(size = 12, family = "ArialMT"),
          axis.text.x  = element_text(size = 12, family = "ArialMT"),
          axis.title.x = element_text(size = 14, family = "ArialMT"),
          legend.title = element_text(size = 14, family = "ArialMT"),
          legend.text = element_text(size = 12, family = "ArialMT"),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 12, family = "ArialMT"))

dev.off()

# Test which SNPs get removed by each problematic loci file
# Load all observed alleles for all patients tested
samples_tested <- fread('./tmp/samples.txt', header = FALSE)$V1

all_loci <- data.table(V1=character(), V2=numeric())
for (this_sample in samples_tested) {
    sample_loci <- list.files(paste0('results/', this_sample, '/snp_input'), full.names = TRUE, pattern = 'loci')
    these_loci <- do.call('rbind', lapply(sample_loci, fread, header = FALSE))
    all_loci <- rbind(all_loci, these_loci)
}

# Remove duplicates
all_loci <- unique(all_loci)
colnames(all_loci) <- c('chr', 'pos')

# Ensure compatible names
all_loci[, chr := gsub('Y','24',gsub('X','23', gsub('chr', '', chr)))]

# Test which SNPs get removed by different probloci files
probloci_removed <- merge(all_loci, filt_probloci_421, all.x = TRUE, by = c('chr', 'pos'))
probloci_removed <- merge(probloci_removed, filt_probloci_bat, all.x = TRUE, by = c('chr', 'pos'))

# Keep only removed SNPs
probloci_removed <- probloci_removed[!(is.na(source.x) & is.na(source.y))]
probloci_removed[, removed_by := ifelse(!is.na(source.x) & !is.na(source.y), 'both',
                                        ifelse(source.x %in% 'tx421', 'tx421', 'battenberg'))]

# Count
nb_probloci_removed <- probloci_removed[, .(nb_removed = .N), by = c('chr', 'removed_by')]

# Get total per caller removing chr Y (24), only present in 421 calls
nb_removed_total <- nb_probloci_removed[chr != 24, .(total_nb = sum(nb_removed)), by = removed_by]

pdf('results/compare_removed_probloci_snps.pdf')
ggplot(nb_probloci_removed, aes(x = chr, y = nb_removed, fill = removed_by)) +
    geom_bar(stat = 'identity', position = position_dodge()) +
    theme_bw() +
    ylab('# SNPs removed') +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size = 14, family = "ArialMT"),
          axis.text.y  = element_text(size = 12, family = "ArialMT"),
          axis.text.x  = element_text(size = 12, family = "ArialMT"),
          axis.title.x = element_text(size = 14, family = "ArialMT"),
          legend.title = element_text(size = 14, family = "ArialMT"),
          legend.text = element_text(size = 12, family = "ArialMT"),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 12, family = "ArialMT"))

ggplot(nb_removed_total, aes(x = removed_by, y = total_nb, fill = removed_by)) +
    geom_bar(stat = 'identity', position = position_dodge()) +
    theme_bw() +
    ylab('# SNPs removed') + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size = 14, family = "ArialMT"),
          axis.text.y  = element_text(size = 12, family = "ArialMT"),
          axis.text.x  = element_text(size = 12, family = "ArialMT"),
          axis.title.x = element_text(size = 14, family = "ArialMT"),
          legend.title = element_text(size = 14, family = "ArialMT"),
          legend.text = element_text(size = 12, family = "ArialMT"),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 12, family = "ArialMT"))

dev.off()