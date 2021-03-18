#!/usr/bin/ python
import sys
import os
import os.path
import csv
from Bio.Seq import Seq
import pandas as pd

########################################################################################################
###                                         scGenotyping.py                                          ###
###         This script extracts the 10x Genomics cell barcode and UMI from PacBio reads             ###
###         derived from sequencing unfragmented cDNA captured by the 10x Genomics protocol.         ###
###                                                                                                  ###
###         This script also assigns a genotype to reads that contain the site of interest.          ###
###                                                                                                  ###
###         To run this script, provide a path to PacBio fastq file, a path to the 10x Genomics      ###
###         cell barcode whitelist, and an output path. The script will output: (1) a text file      ###
###         containing each read's cell barcode, UMI, and transcript sequence, (2) a narrowed        ###
###         down version of the first file that only contains unique reads and valid cell barcodes,  ###
###         and (3) a set of files containing the genotype annotations for each of those reads.      ###
###                                                                                                  ###
###                                     Bedrosian Lab, July 2020                                     ###
########################################################################################################

def processPacBio(pacbio_fastq, whitelist, out_path):
    print("Processing " + pacbio_fastq + " to " + out_path + " using specified whitelist " + whitelist)

    # Beginning with original PacBio fastq file, subset the reads containing a partial 10x Genomics Read 1 sequence,
    # and output a text file containing the cell barcode, umi, and transcript sequence::
    with open(pacbio_fastq) as fastq, open(out_path + "fulloutput.txt", 'wb') as outfile:

            primer_fwd = "CGACGCTCTTCCGATCT"
            primer_rev = "AGATCGGAAGAGCGTCG"
            line_count = 0

            for line in fastq:
                line = line.strip()
                index_fwd = line.find(primer_fwd)
                index_rev = line.find(primer_rev)
                if index_fwd >= 0:
                    barcode = line[index_fwd + len(primer_fwd):index_fwd + len(primer_fwd) + 16]
                    umi = line[index_fwd + len(primer_fwd) + 16: index_fwd + len(primer_fwd) + 28]
                    transcript = line[index_fwd + len(primer_fwd) + 28: -1]
                    line_count += 1
                    outfile.write(barcode + "\t" + umi + "\t" + transcript + "\n")
                elif index_rev >= 0:
                    seq = str(Seq(line).reverse_complement())
                    index_fwd_revcomp = seq.find(primer_fwd)
                    barcode = seq[index_fwd_revcomp + len(primer_fwd):index_fwd_revcomp + len(primer_fwd) + 16]
                    umi = seq[index_fwd_revcomp + len(primer_fwd) + 16: index_fwd_revcomp + len(primer_fwd) + 28]
                    transcript = seq[index_fwd_revcomp + len(primer_fwd) + 28: -1]
                    line_count += 1
                    outfile.write(barcode + "\t" + umi + "\t" + transcript + "\n")
            print("Number of reads containing partial 10x read 1 sequence: " + str(line_count))

    # Cross-reference observed cell barcodes against 10x Genomics whitelist of known barcodes to output reads with legitimate
    # cell barcodes. Remove UMI duplicates.
    with open(out_path + "fulloutput.txt", 'rb') as seenbarcodes, open(whitelist, 'rb') as knownbarcodes:
        with open(out_path + "goodbarcodes.txt", 'wb') as goodbarcodes:

                df1 = pd.read_csv(seenbarcodes, sep = '\t', names = ['barcode', 'umi', 'transcript'])
                df2 = pd.read_csv(knownbarcodes, sep = '\t', names = ['barcode'])
                df3 = pd.merge(df1, df2, on = 'barcode', how = 'inner')
                df4 = df3.drop_duplicates('umi')
                print("Number of unique reads with legitimate 10x cell barcodes: " + str(len(df4.index)))
                df4.to_csv(goodbarcodes, sep='\t', index=False)

    # Genotype reads for first or second PTEN indel.
    with open(out_path + "goodbarcodes.txt", 'rb') as file1:
        with open (out_path + "pten1_geno.txt", 'wb') as outfile:

                pten1_fwd = "TCCTTTTGAAGACCATAACC"
                pten1_rev = "GGTTATGGTCTTCAAAAGGA"
                line_count = 0
                line_count_wt = 0
                line_count_mut = 0

                for line in file1:
                    line = line.strip()
                    index_fwd = line.find(pten1_fwd)
                    index_rev = line.find(pten1_rev)
                    if index_fwd >= 0:
                        genotype = line[index_fwd-3:index_fwd]
                        wt = "ATA"
                        mut = "TCT"
                        if wt in genotype:
                            status = "WT"
                            line_count_wt += 1
                        elif mut in genotype:
                            status = "MUT"
                            line_count_mut += 1
                        else:
                            status = "Not defined"
                        outfile.write(line[:-1] + "\t" + genotype + "\t" + status + "\n")
                        line_count += 1
                    elif index_rev >= 0:
                            genotype = line[index_rev + len(pten1_rev):index_rev + len(pten1_rev) + 3]
                            wt = "TAT"
                            mut = "AGA"
                            if wt in genotype:
                                status = "WT"
                                line_count_wt += 1
                            elif mut in genotype:
                                status = "MUT"
                                line_count_mut += 1
                            else:
                                status = "Not defined"
                            outfile.write(line[:-1] + "\t" + genotype + "\t" + status + "\n")
                            line_count += 1
                print("PTEN1 genotyped reads: " + str(line_count))
                print("PTEN1 wildtype reads: " + str(line_count_wt))
                print("PTEN1 mutant reads: " + str(line_count_mut))

    with open(out_path + "goodbarcodes.txt", 'rb') as file1:
        with open (out_path + "pten2_geno.txt", 'wb') as outfile:

                pten2_fwd = "TTCTGTAACACCAGATGTTA"
                pten2_rev = "TAACATCTGGTGTTACAGAA"
                line_count = 0
                line_count_wt = 0
                line_count_mut = 0

                for line in file1:
                    line = line.strip()
                    index_fwd = line.find(pten2_fwd)
                    index_rev = line.find(pten2_rev)
                    if index_fwd >= 0:
                        genotype = line[index_fwd + len(pten2_fwd):index_fwd + len(pten2_fwd) + 5]
                        wt = "GTGAC"
                        mut = "GTGTG"
                        if wt in genotype:
                            status = "WT"
                            line_count_wt += 1
                        elif mut in genotype:
                            status = "MUT"
                            line_count_mut += 1
                        else:
                            status = "Not defined"
                        outfile.write(line[:-1] + "\t" + genotype + "\t" + status + "\t" + "\n")
                        line_count += 1
                    elif index_rev >= 0:
                        seq = str(Seq(line).reverse_complement())
                        index_fwd_revcomp = seq.find(pten2_fwd)
                        if index_fwd_revcomp >= 0:
                            genotype = seq[index_fwd_revcomp + len(pten2_fwd):index_fwd_revcomp + len(pten2_fwd) + 5]
                            wt = "GTGAC"
                            mut = "GTGTG"
                            if wt in genotype:
                                status = "WT"
                                line_count_wt += 1
                            elif mut in genotype:
                                status = "MUT"
                                line_count_mut += 1
                            else:
                                status = "Not defined"
                            outfile.write(line[:-1] + "\t" + genotype + "\t" + status + "\t" + "\n")
                            line_count += 1
                print("PTEN2 genotyped reads: " + str(line_count))
                print("PTEN2 wildtype reads: " + str(line_count_wt))
                print("PTEN2 mutant reads: " + str(line_count_mut))

n_arguments = len(sys.argv)

if n_arguments > 3:
    infile1 = sys.argv[1]
    infile2 = sys.argv[2]
    outbase = sys.argv[3]
    files = [infile1, infile2]
    exists = [f for f in files if os.path.isfile(f)];
    non_exist = list(set(exists) ^ set(files))

    print("existing: %s" % exists)
    print("non existing: %s" % non_exist)

    if exists:
        processPacBio(infile1, infile2, outbase)
    else:
        print("ERROR: Missing file")
        print("USAGE: python " + sys.argv[0] + " PacBioFastq 10xWhitelist OutPath")
        sys.exit(-1)
else:
    print("USAGE: python " + sys.argv[0] + " PacBioFastq 10xWhitelist OutPath")
    print("USAGE: Please provide a path to PacBio fastq file, path to 10x Genomics cell barcode whitelist, and path to output directory")
    sys.exit(-1)
