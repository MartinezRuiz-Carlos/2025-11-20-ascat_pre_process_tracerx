library(data.table)
library(ggplot2)
library(GenomicRanges)
library(gtools)
library(argparse)

# Script to test differences between previous iterations of logR and BAF on varscan compared with the new methods

parser <- ArgumentParser(description = "Compare ASCAT HTS BAF and logR with previous runs based on Varscan")
parser$add_argument("--varscan_obj_path", help = "Path to Varscan output object")
parser$add_argument("--sample_name_hash", help = "Sample name followed by sequencing hash")
parser$add_argument("--ascat_baf",        help = "Path to ASCAT HTS BAF")
parser$add_argument("--ascat_logr",       help = "Path to ASCAT HTS logR")
parser$add_argument("--results_dir",      help = "Directory to store results")

args <-parser$parse_args()

# Load varscan BAF and logr and filter -----------------------------------------------------
varscan_object <- readRDS(args$varscan_obj_path)

# Get chromosomes in the right format
varscan_pos <- varscan_object$SNPpos
varscan_pos$chrom <- paste0('chr', varscan_pos$chrom)

# Bind together
varscan_baf <- as.data.table(cbind(varscan_pos, varscan_object$Tumor_BAF))
varscan_logr <- as.data.table(cbind(varscan_pos, varscan_object$Tumor_LogR))

# Rename colnames and filter out NAs
colnames(varscan_baf) <- c('chrom', 'pos', 'baf')
varscan_baf <- varscan_baf[!is.na(baf)]
varscan_baf[, source := 'varscan']

colnames(varscan_logr) <- c('chrom', 'pos', 'logr')
varscan_logr <- varscan_logr[!is.na(logr)]
varscan_logr[, source := 'varscan']

# Load ASCAT generated BAF and logR files---------------------------------------------------------
ascat_baf <- fread(args$ascat_baf)
ascat_logr <- fread(args$ascat_logr)

# Keep columns of interest and ensure compatibility with varscan liftover
cols_keep_ascat <- c('Chromosome', 'Position', args$sample_name_hash)
ascat_baf_filt <- ascat_baf[, ..cols_keep_ascat]
setnames(ascat_baf_filt, cols_keep_ascat, c('chrom', 'pos', 'baf'))
ascat_baf_filt[, chrom := paste0('chr', chrom)]
ascat_baf_filt[, source := 'ascat']

ascat_logr_filt <- ascat_logr[, ..cols_keep_ascat]
setnames(ascat_logr_filt, cols_keep_ascat, c('chrom', 'pos', 'logr'))
ascat_logr_filt[, chrom := paste0('chr', chrom)]
ascat_logr_filt[, source := 'ascat']

# Merge and plot----------------------------------------------------------------------------------------
baf_compare <- rbind(varscan_baf, ascat_baf_filt)
logr_compare <- rbind(varscan_logr, ascat_logr_filt)

pdf(paste0(args$results_dir, '/plot_comparisons/logr_baf_diffs.pdf'), height = 6, width = 22)
# Sort chromosomes
baf_compare[, chr_factor := factor(chrom, levels = mixedsort(unique(chrom)))]
logr_compare[, chr_factor := factor(chrom, levels = mixedsort(unique(chrom)))]
ggplot(baf_compare, aes(x = pos, y = baf, colour = source)) +
    geom_point(size = 0.1, alpha = 0.1) +
    facet_grid(. ~ chr_factor, scales = 'free_x') +
    geom_vline(data = baf_compare[chr_factor == 'chr6'], aes(xintercept = 28477797), colour = 'blue', size = 0.1) + # Mark MHC location
    geom_vline(data = baf_compare[chr_factor == 'chr6'], aes(xintercept = 33448354), colour = 'blue', size = 0.1) +
    theme_bw() +
    ylab('BAF') + xlab('Position') +
    guides(color = guide_legend(override.aes = list(size = 1, alpha = 1))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size = 14, family = "ArialMT"),
          axis.text.y  = element_text(size = 12, family = "ArialMT"),
          axis.ticks.x =element_blank(),
          axis.text.x  = element_blank(),
          axis.title.x = element_text(size = 14, family = "ArialMT"),
          legend.title = element_text(size = 14, family = "ArialMT"),
          legend.text = element_text(size = 12, family = "ArialMT"),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 12, family = "ArialMT"))

ggplot(logr_compare, aes(x = pos, y = logr, colour = source)) +
    geom_point(size = 0.1, alpha = 0.1) +
    facet_grid(. ~ chr_factor, scales = 'free_x') +
    geom_vline(data = logr_compare[chr_factor == 'chr6'], aes(xintercept = 28477797), colour = 'blue', size = 0.1) + # Mark MHC location
    geom_vline(data = logr_compare[chr_factor == 'chr6'], aes(xintercept = 33448354), colour = 'blue', size = 0.1) +
    theme_bw() +
    ylab('logR') + xlab('Position') +
    guides(color = guide_legend(override.aes = list(size = 1, alpha = 1))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size = 14, family = "ArialMT"),
          axis.text.y  = element_text(size = 12, family = "ArialMT"),
          axis.ticks.x =element_blank(),
          axis.text.x  = element_blank(),
          axis.title.x = element_text(size = 14, family = "ArialMT"),
          legend.title = element_text(size = 14, family = "ArialMT"),
          legend.text = element_text(size = 12, family = "ArialMT"),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 12, family = "ArialMT"))
dev.off()

# Average BAF and logR in 1kb windows for a "fairer" comparison
window_len <- 10000
baf_compare_avg <- baf_compare[, .(baf_avg = mean(baf, na.rm = TRUE)), 
                               by = .(chrom, window_start = floor(pos / window_len) * window_len, source)]

logr_compare_avg <- logr_compare[, .(logr_avg = mean(logr, na.rm = TRUE)), 
                               by = .(chrom, window_start = floor(pos / window_len) * window_len, source)]

pdf(paste0(args$results_dir, '/plot_comparisons/logr_baf_diffs_avg.pdf'), height = 6, width = 22)
baf_compare_avg[, chr_factor := factor(chrom, levels = mixedsort(unique(chrom)))]
logr_compare_avg[, chr_factor := factor(chrom, levels = mixedsort(unique(chrom)))]

ggplot(baf_compare_avg, aes(x = window_start, y = baf_avg, colour = source)) +
    geom_point(size = 0.1, alpha = 0.3) +
    facet_grid(. ~ chr_factor, scales = 'free_x') +
    theme_bw() +
    ylab('BAF') + xlab('Position') +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size = 14, family = "ArialMT"),
          axis.text.y  = element_text(size = 12, family = "ArialMT"),
          axis.ticks.x =element_blank(),
          axis.text.x  = element_blank(),
          axis.title.x = element_text(size = 14, family = "ArialMT"),
          legend.title = element_text(size = 14, family = "ArialMT"),
          legend.text = element_text(size = 12, family = "ArialMT"),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 12, family = "ArialMT"))

ggplot(logr_compare_avg, aes(x = window_start, y = logr_avg, colour = source)) +
    geom_point(size = 0.1, alpha = 0.3) +
    facet_grid(. ~ chr_factor, scales = 'free_x') +
    theme_bw() +
    ylab('logR') + xlab('Position') +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size = 14, family = "ArialMT"),
          axis.text.y  = element_text(size = 12, family = "ArialMT"),
          axis.ticks.x =element_blank(),
          axis.text.x  = element_blank(),
          axis.title.x = element_text(size = 14, family = "ArialMT"),
          legend.title = element_text(size = 14, family = "ArialMT"),
          legend.text = element_text(size = 12, family = "ArialMT"),
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 12, family = "ArialMT"))
dev.off()

# Correlate BAF and logr. There is good overlap with BAF positions, but logR needs to be padded (original varscan is at 100bp intervals, so add 100bp to the right of start and overlap)
baf_corr_dt <- merge(varscan_baf, ascat_baf_filt, by = c('chrom', 'pos'))
setnames(baf_corr_dt, c('baf.x', 'baf.y'), c('baf_varscan', 'baf_ascat'))

varscan_logr_gr <- GRanges(
  seqnames = varscan_logr$chrom,
  ranges = IRanges(start = varscan_logr$pos, end = varscan_logr$pos + 100),
  logr_varscan = varscan_logr$logr
)

ascat_logr_gr <- GRanges(
  seqnames = ascat_logr_filt$chrom,
  ranges = IRanges(start = ascat_logr_filt$pos, end = ascat_logr_filt$pos),
  logr_ascat = ascat_logr_filt$logr
)

logr_corr_dt <- unique(as.data.table(mergeByOverlaps(varscan_logr_gr, ascat_logr_gr)))

# Plot
pdf(paste0(args$results_dir, '/plot_comparisons/logr_baf_corrs.pdf'), width = 8, height = 6)
baf_corr_dt[, chr_type := ifelse(chrom == 'chrX', 'X', 'Autosomes')]
ggplot(baf_corr_dt, aes(x = baf_varscan, y = baf_ascat)) +
    geom_bin2d(bins = 1000) +
    geom_abline(colour = 'red') +
    geom_abline(slope = -1, intercept = 1, colour = 'red') +
    theme_bw() +
    xlim(0, 1) + ylim(0, 1) +
    facet_grid(. ~ chr_type) +
    ylab('BAF Varscan') + xlab('BAF ASCAT') +
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

logr_corr_dt[, chr_type := ifelse(varscan_logr_gr.seqnames == 'chrX', 'X', 'Autosomes')]
ggplot(logr_corr_dt, aes(x = logr_varscan, y = logr_ascat)) +
    geom_bin2d(bins = 1000) + 
    theme_bw() +
    xlim(-3.7, 3.7) + ylim(-3.7, 3.7) +
    geom_abline(colour = 'red') +
    facet_grid(. ~ chr_type) + 
    ylab('logR Varscan') + xlab('logR ASCAT') +
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