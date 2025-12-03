import pysam
import pandas as pd
import os
import re
import argparse

# Script to combine 1000G SNP file with actual SNP calls per patient and turn them into the ascat_HTS input format

# Functions
def skip_variant(record, min_dp, min_ad, caller):
    """
    Function to test if a variant needs to be skipped, will skipped variants if
        - Depth < min_dp
        - Variant count < min_ad
        - The variant is homozygote
        - The variant is not PASS
        - The variant is an indel (reference/alt alleles are longer than 1 nucleotide)

    Parameters:
        - record: pysam.VariantRecord. SNP record
        - min_dp: int. Minimum depth to filter by
        - min_ad: int. Minimum alternative allele to filter by
        - caller: str. Tool used for SNP calling, one of deepvariant or platypus

    Returns:
        Boolean as to whether or not the SNP should be skipped (filtered)

    """
    # Get the depth (DP), genotype (GT) and variant count (AD) for each variant, assume only one sample
    this_sample = record.samples[0]
    if caller == 'deepvariant':
        depth = this_sample['DP'][0]
        # Keep only last value (only bi-allelic variants are going to be kept, so we are confident that will always be the right value)
        variant_count = this_sample['AD'][1]
    elif caller == 'platypus': 
        depth = this_sample['NR'][0]
        variant_count = this_sample['NV'][0]
    else:
        raise Exception('caller needs to be one of "platypus" or "deepvariant"')
    gt = this_sample['GT']
    #Filter
    filt = list(record.filter)[0]
    # Ref alt to check if variant is bi-allelic
    ref = record.ref
    alt = record.alts[0]

    # Get all conditions
    single_nts = ['A', 'T', 'C', 'G']
    # Depth les than threshold, varcount less than threshold, homozygous variant, not PASS, not single nt reference (indel), not bi-allelic, not single nt alt (indel)
    skip_var = (depth < min_dp) or (variant_count < min_ad) or (gt[0] == gt[1]) or (filt != 'PASS') or (ref not in single_nts) or (len(record.alts) > 1) or (alt not in single_nts)
    
    return(skip_var)

def format_nts(nt_list):
    """
    Function format ref and alt nucleotides to ASCAT format A=1, C=2, G=3 and T=4

    Parameters:
        - nt_list: list. List of nucleotides to re-format
    
    Returns:
        - List of reformatted nucleotides for ASCAT

    """
    nt_convert = {'A':1, 'C':2, 'G':3, 'T':4}
    nt_out = [nt_convert[nt] for nt in nt_list]
    return(nt_out)


def main():
    parser = argparse.ArgumentParser(description='Script to combine 1000G SNPs with het patient-specific SNPs and parse into ASCAT format')
    parser.add_argument('-i', '--input_vcf', help = 'Input germline SNP VCF - Currently expects a Platypus VCF', required = True)
    parser.add_argument('-c', '--caller', help = 'Variant caller used to generate the VCF file, one of "deepvariant" or "platypus"', required = True)    
    parser.add_argument('-l', '--loci_input', help = 'Path to directory with 1000G alleles formatted for ASCAT - Assumes the chr name is on the file before .txt', required = True)
    parser.add_argument('-d', '--min_dp', help = 'Minimum depth to filter SNPs by', default = 30)
    parser.add_argument('-a', '--min_ad', help = 'Minimum alt allele depth per SNP', default = 10)
    parser.add_argument('-o', '--out_dir', default='.', help = 'Output directory')

    args = parser.parse_args()

    # Get patient SNP VCF
    vcf = pysam.VariantFile(args.input_vcf, 'r')

    # Extract position, ref/alt, ONLY if the variant passes quality and depth filters
    all_snps = list()
    for record in vcf.fetch():

        if (skip_variant(record, args.min_dp, args.min_ad, args.caller)):
            continue
        else:
            snp_dict = {'chr':record.chrom, 'position':record.pos, 'ref':record.ref, 'alt':record.alts[0]}
            all_snps.append(snp_dict)

    all_snps_df = pd.DataFrame(all_snps)

    # Load 1000G loci files as a single data frame
    all_allele_files = [os.path.join(args.loci_input, file) for file in os.listdir(args.loci_input)]

    all_allele_1000g_list = []
    for allele_file in all_allele_files:
        chr_allele = pd.read_csv(allele_file, sep = '\t')
        # Get chromosome name from file name
        chr_name = re.sub(r'(.+_)(chr.+)(\.txt)', r'\2', allele_file)
        chr_allele['chr'] = chr_name

        all_allele_1000g_list.append(chr_allele)

    all_allele_1000g_df = pd.concat(all_allele_1000g_list, ignore_index = True)
    all_allele_1000g_df['snp_label'] = all_allele_1000g_df['chr'] + ':' + all_allele_1000g_df['position'].astype(str)

    # Remove the positions that are already present in the patient VCF
    patient_snps = (all_snps_df['chr'] + ':' + all_snps_df['position'].astype(str)).tolist()
    all_allele_1000g_df = all_allele_1000g_df[~all_allele_1000g_df['snp_label'].isin(patient_snps)]

    # Convert nucleotides into the format ASCAT expects for the patient SNPs and merge
    all_snps_df['a0'] = format_nts(all_snps_df['ref'])
    all_snps_df['a1'] = format_nts(all_snps_df['alt'])

    # Get the columns of interest and merge
    all_allele_1000g_df = all_allele_1000g_df[['chr', 'position', 'a0', 'a1']]
    joint_snps_df = all_snps_df[['chr', 'position', 'a0', 'a1']].copy()

    joint_snps_df = pd.concat([joint_snps_df, all_allele_1000g_df])

    # Sort
    joint_snps_df.sort_values(by = ['position'],
                            inplace = True)

    # Write allele and loci files to output
    for this_chr in joint_snps_df['chr'].drop_duplicates().tolist():
        out_chr_df = joint_snps_df[joint_snps_df['chr'] == this_chr]
        
        out_chr_alleles = out_chr_df[['position', 'a0', 'a1']]
        out_chr_loci = out_chr_df[['chr', 'position']]

        out_file_alleles = os.path.join(args.out_dir, 'alleles_' + this_chr + '.txt')
        out_file_loci = os.path.join(args.out_dir, 'loci_' + this_chr + '.txt')

        out_chr_alleles.to_csv(out_file_alleles, sep = '\t',
                            index = False, header = True,
                            na_rep = '', quoting = 3)

        out_chr_loci.to_csv(out_file_loci, sep = '\t',
                            index = False, header = False,
                            na_rep = '', quoting = 3)
        
if __name__ == "__main__":
        main()
