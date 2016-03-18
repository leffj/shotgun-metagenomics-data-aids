#!/usr/bin/env python

__author__ = "Jonathan Leff, Hannah Holland-Moritz"
__email__ = "jonathan.leff@colorado.edu"
__version__ = "0.0.1"

"""
Create an OTU table from a set of shotgun metagenomes using metaxa and usearch.
Requires: cutadapt, USEARCH v7, metaxa2
metaxa2 executable must be named 'metaxa2'
"""

import argparse
import gzip
import subprocess
import os
from collections import defaultdict
from collections import Counter


def main():
    parser = argparse.ArgumentParser(description=\
        'Create otu table from shotgun files using metaxa2 and usearch. Make sure \
        metaxa2 is installed properly before running.')
    req = parser.add_argument_group('required arguments')
    req.add_argument('-i', '--input_fp', required=True,
        type=str, help='A tab-delimited file with sample IDs in the first \
        column, corresponding filepaths for the forward reads and reverse \
        reads in the second and third columns, respectively. No header or \
        header must start with "#".')
    req.add_argument('-d', '--database_fp', required=True,
        help='The database filepath. So far, only tested with greengenes \
        clustered at 97%% similarity. Format is fasta.')
    parser.add_argument('-o', '--output_dir',
        default='OTU_table_from_shotgun_data_files',
        help='The output directory. Will create if does not exist. In addition \
        to the OTU table, intermediate files will be produced in this dir.')
    parser.add_argument('-a', '--adapter_sequence', default='CTGTCTCTTATA',
        help='The adapter sequence that will be trimmed from both forward \
        and reverse paired-end reads. The default is the standard Nextera \
        adapter.')
    parser.add_argument('-p', '--processors', default=1,
        help='The number of processors to use when running on a multicore \
        computer.')
    parser.add_argument('-t', '--taxonomy_fp',
        help='The filepath for the taxonomy file that links the database \
        identifiers with taxonomy strings. If provided, taxonomy strings \
        will be added to the output OTU table.')
    parser.add_argument('-g', '--taxonomic_group', default='bacteria',
        help='Use sequences extracted for either: "bacteria", "archaea", or \
        "eukaryota". Default = "bacteria"')
    parser.add_argument('-c', '--combined_reads', default=False,
        help='Set this parameter to "True" if paired-end reads are \
        interleaved (i.e. they alternate in one file). Uses "pefcon" to \
        split the reads into separate files.')
    parser.add_argument('-m', '--metaxa_sequences_dir',
        help='This argument can be used to specify the metaxa output directory \
        of a previous run of this script so you do not have to start from the \
        beginning if investigating another taxonomic group or the run crashes.')
    # parser.add_argument('-r', '--merge_reads', default=True,
    #     help='Set this parameter to "False" to disable merging of paired-end \
    #     reads. Default = True.') # to do this, need to figure out an alternate
    # way to map reads to db, since the 'N's that metaxa creates are problematic
    # for usearch


    args = parser.parse_args()

    # create output directory if doesn't exist
    out_dir = args.output_dir
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    else:
        raise RuntimeError('Output directory already exists. Remove or specify a \
            new directory.')
    out_dir = out_dir + "/"


    # get sample list and associated input files. Deal with interleaved samples
    # if needed
    if args.combined_reads:
        samples, R_fps = parse_input_file(args.input_fp, interleaved=True)
        R1_fps = {}
        R2_fps = {}
        for sample in samples:
            print "Spliting interleaved reads for " + sample + "..."
            R1_fps[sample], R2_fps[sample] = \
                split_interleaved_reads(sample, R_fps[sample], out_dir)
    else:
        samples, R1_fps, R2_fps = parse_input_file(args.input_fp)

    # start with metaxa seqs
    otu_counts = {}
    if args.metaxa_sequences_dir:
        for sample in samples:
            metaxa_extract_fp = "%s/%s_metaxa.%s.fasta" % \
                (args.metaxa_sequences_dir, sample, args.taxonomic_group)
            print metaxa_extract_fp
            print "Mapping %s sequences to database for %s..." % \
                (args.taxonomic_group, sample)
            readmap_fp = map_seqs_to_db(sample, metaxa_extract_fp, args.database_fp,
                                     args.processors, out_dir)
            print "Generating OTU counts for " + sample + "..."
            ssu_seq_count, ssu_hits = return_number_SSU_seqs_hits(readmap_fp)
            otu_counts[sample] = convert_readmap_to_OTU_counts(sample,
                                                               readmap_fp,
                                                               out_dir)
        print "Writting OTU table..."
        merge_otu_counts_and_write(otu_counts, out_dir + 'otu_table.txt')

        # optionally, add taxonomy strings
        if args.taxonomy_fp:
            print "Adding taxonomic classifications..."
            add_taxonomy_to_otu_table(out_dir + 'otu_table.txt', args.taxonomy_fp)
        else:
            print "No taxonomy file provided."

    # or start from beginning
    else:
        # prep summary file
        summary_stats_out = open(out_dir + "summary_stats.txt", "w")
        summary_stats_out.write("#Sample_ID\tRaw_seq_count\tMerged_seq_count\t\
            Mean_merged_seq_length\tSSU_seq_count\tSSU_hits_database\tYield\n")

        # perform analyses for each sample and record summary stats along the way
        for sample in samples:
            print "Removing adapters from " + sample + "..."
            R1_trim_fp, R2_trim_fp = remove_adapters(sample, R1_fps[sample],
                                                     R2_fps[sample],
                                                     args.adapter_sequence,
                                                     out_dir)
            print "Merging paired reads from " + sample + "..."
            merged_fp = merge_paired_reads(sample, R1_trim_fp, R2_trim_fp, out_dir)
            print "Getting sequence stats from " + sample + "..."
            merged_seq_count, merged_seq_length = get_seq_stats(merged_fp)
            print "Finding SSU rRNA sequences from " + sample + "..."
            metaxa_extract_fp = run_metaxa(sample, merged_fp, args.taxonomic_group,
                                           args.processors, out_dir)
            print "Mapping %s sequences to database for %s..." % \
                    (args.taxonomic_group, sample)
            readmap_fp = map_seqs_to_db(sample, metaxa_extract_fp, args.database_fp,
                                     args.processors, out_dir)
            print "Generating OTU counts for " + sample + "..."
            ssu_seq_count, ssu_hits = return_number_SSU_seqs_hits(readmap_fp)
            otu_counts[sample] = convert_readmap_to_OTU_counts(sample,
                                                               readmap_fp,
                                                               out_dir)
            # Write summary statistics to file
            if R1_fps[sample].endswith('.gz'):
                seqs_opened = gzip.open(R1_fps[sample], 'rb')
            else:
                seqs_opened = open(R1_fps[sample], 'U')
            summary_stats_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                                    %(sample, get_seq_count(
                                    seqs_opened),
                                    merged_seq_count, merged_seq_length,
                                    ssu_seq_count, ssu_hits,
                                    (float(ssu_hits) / float(merged_seq_count))))

        # merge OTU sequence count data across samples and write to file
        # otu_counts = {'test1': {'289883': 1, '4469420': 1, '4402605': 1, '4425441': 2,
        #               '851811': 1, '4470260': 1},
        #               'test2': {'1138744': 1, '4425441': 5}}
        print "Writting OTU table..."
        merge_otu_counts_and_write(otu_counts, out_dir + 'otu_table.txt')

        # optionally, add taxonomy strings
        if args.taxonomy_fp:
            print "Adding taxonomic classifications..."
            add_taxonomy_to_otu_table(out_dir + 'otu_table.txt', args.taxonomy_fp)
        else:
            print "No taxonomy file provided."



def parse_input_file(input_fp, interleaved=False):
    input = open(input_fp, 'U')
    sampleIDs = []
    if interleaved:
        sample_fps = {}
        for line in input:
            fields = line.strip().split('\t')
            if fields[0].startswith('#'):
                continue
            sampleIDs.append(fields[0])
            # check if filepath exists
            if not os.path.exists(fields[1]):
                raise RuntimeError("Filepath does not exist: %s" %(fields[1]))
            sample_fps[fields[0]] = fields[1]
        return sampleIDs, sample_fps
    else:
        sample_fps_R1s = {}
        sample_fps_R2s = {}
        for line in input:
            fields = line.strip().split('\t')
            if fields[0].startswith('#'):
                continue
            sampleIDs.append(fields[0])
            # check if filepath exists
            if not os.path.exists(fields[1]):
                raise RuntimeError("Filepath does not exist: %s" %(fields[1]))
            if not os.path.exists(fields[2]):
                raise RuntimeError("Filepath does not exist: %s" %(fields[2]))
            sample_fps_R1s[fields[0]] = fields[1]
            sample_fps_R2s[fields[0]] = fields[2]
        return sampleIDs, sample_fps_R1s, sample_fps_R2s

def remove_adapters(sampleID, R1_in_fp, R2_in_fp, adapter, out_dir):
    if not os.path.exists(out_dir + "adapter_removal"):
        os.makedirs(out_dir + "adapter_removal")
    R1_out_fp = out_dir + "adapter_removal/" + sampleID + "_R1_trim.fq"
    R2_out_fp = out_dir + "adapter_removal/" + sampleID + "_R2_trim.fq"
    info_out_fp = out_dir + "adapter_removal/" + sampleID + "_cutadapt_info.txt"
    out_summary_fp = out_dir + "adapter_removal/" + sampleID + "_cutadapt_summary.txt"
    out_summary = open(out_summary_fp, "w")
    subprocess.call(["cutadapt", "-a", adapter, "-A", adapter, "-o",
                    R1_out_fp, "-p", R2_out_fp, R1_in_fp, R2_in_fp],
                    stdout = out_summary)
    out_summary.close()
    return R1_out_fp, R2_out_fp

def merge_paired_reads(sampleID, R1_fp, R2_fp, out_dir):
    if not os.path.exists(out_dir + "merged_paired_reads"):
        os.makedirs(out_dir + "merged_paired_reads")
    merged_fp = out_dir + "merged_paired_reads/" + sampleID + "_merged.fq"
    FNULL = open(os.devnull, 'w')
    merge_log = open(out_dir + "merged_paired_reads/" + sampleID + "_merge_log.txt", "w")
    subprocess.call(["usearch7", "-fastq_mergepairs", R1_fp, "-reverse",
                    R2_fp, "-fastqout", merged_fp, "-fastq_truncqual", "3",
                    "-fastq_maxdiffs", "5", "-fastq_minovlen", "10",
                    "-fastq_minmergelen", "50"], stdout = FNULL,
                    stderr = merge_log)
    return merged_fp

def get_seq_stats(input_fq_fp):
    input_fq = open(input_fq_fp)
    seq_count = 0
    seq_lengths = []
    for header, sequence, quality in basic_fastq_parser(input_fq):
        seq_count += 1
        seq_lengths.append(len(sequence))
    mean_seq_len = round(sum(seq_lengths)/float(len(seq_lengths)), 2)
    return seq_count, mean_seq_len

def basic_fastq_parser(in_f):
    lineno, head, seq, qual = 0, "", "", ""
    for l in in_f:
        lineno += 1
        if lineno % 4 == 1:
            head = l.strip()
        elif lineno % 4 == 2:
            seq = l.strip()
        elif lineno % 4 == 0:
            qual = l.strip()
            yield head, seq, qual

def get_seq_count(in_f):
    count = 0
    for head, seq, qual in basic_fastq_parser(in_f):
        count += 1
    return count

def run_metaxa(sampleID, input_fq_fp, domain, procs, out_dir):
    if not os.path.exists(out_dir + "metaxa_analysis"):
        os.makedirs(out_dir + "metaxa_analysis")
    output_pre = out_dir + "metaxa_analysis/" + sampleID + "_metaxa"
    metaxa_log_fp = open(out_dir + "metaxa_analysis/" + sampleID + "_metaxa_stdout_log", "w")
    subprocess.call(["metaxa2", "-i", input_fq_fp, "-plus", "T", "--cpu",
    str(procs), "-o", output_pre], stderr = metaxa_log_fp)
    if os.path.exists("error.log"):
        subprocess.call(["mv", "error.log", out_dir + "metaxa_analysis/" + \
                         sampleID + "error.log"])
    return output_pre + "." + domain + ".fasta"

def map_seqs_to_db(sampleID, input_fp, db_fp, procs, out_dir):
    if not os.path.exists(out_dir + "reads_mapped_to_db"):
        os.makedirs(out_dir + "reads_mapped_to_db")
    out_fp = out_dir + "reads_mapped_to_db/" + sampleID + "_readmap.uc"
    FNULL = open(os.devnull, 'w')
    map_log = open(out_dir + "reads_mapped_to_db/" + sampleID + "_readmap_log.txt", "w")
    subprocess.call(["usearch7", "-usearch_global", input_fp, "-db", db_fp,
                    "-id", "0.97", "-strand", "both", "-uc", out_fp,
                    "-maxrejects", "0", "-threads", str(procs)], stdout = FNULL,
                    stderr = map_log)
    return out_fp

def return_number_SSU_seqs_hits(readmap_fp):
    number_seqs = 0
    number_hits = 0
    for line in open(readmap_fp):
        number_seqs += 1
        if line.split('\t')[0] == 'H':
            number_hits += 1
    return number_seqs, number_hits

def convert_readmap_to_OTU_counts(sampleID, readmap_fp, out_dir):
    otu_counts = {}
    readmap = open(readmap_fp, 'r')
    for line in readmap:
        line = line.rstrip()
        line_split = line.split('\t')
        if line_split[0] == 'H':
            if line_split[9] in otu_counts:
                otu_counts[line_split[9]] += 1
            else:
                otu_counts[line_split[9]] = 1
    return otu_counts

def merge_otu_counts_and_write(otu_counts_dict, out_fp):
    table = defaultdict(Counter)
    samples = []
    for sample in otu_counts_dict:
        samples.append(sample)
        for otu in otu_counts_dict[sample]:
            table[otu][sample] += otu_counts_dict[sample][otu]
    # write table
    output = open(out_fp, "w")
    line = "#OTU_ID\t%s\n"%('\t'.join(samples))
    output.write(line)
    for otu in table.keys():
        counts = [0] * len(samples)
        for i, sample in enumerate(samples):
            counts[i] = table[otu][sample]
        line = "%s\t%s\n" %(otu, '\t'.join(map(str, counts)))
        output.write(line)
    output.close()

def add_taxonomy_to_otu_table(otu_table_fp, taxonomy_fp):
    otu_table = open(otu_table_fp, 'U')
    taxonomy = open(taxonomy_fp, 'U')
    taxonomy_dict = {}
    table_out = []
    for line in taxonomy:
        line_split = line.strip().split('\t')
        taxonomy_dict[line_split[0]] = line_split[1]
    for line in otu_table:
        line_split = line.strip().split('\t')
        if line_split[0].startswith('#'):
            table_out.append(line.strip() + "\ttaxonomy\n")
        else:
            otu_tax = taxonomy_dict[line_split[0]]
            table_out.append("%s\t%s\n" %(line.strip(), otu_tax))
    otu_table.close()
    output = open(otu_table_fp, 'w')
    output.write("".join(table_out))
    output.close()

def split_interleaved_reads(sampleID, in_fp, out_dir):
    # pefcon -1 /data/shared/2013_06_NutNet-metagenomic/raw_data/6497.7.45694.AATAGG.fastq.gz -o 6497.7.45694.AATAGG -z --split
    if not os.path.exists(out_dir + "split_paired_reads"):
        os.makedirs(out_dir + "split_paired_reads")
    out_pre = out_dir + "split_paired_reads/" + sampleID
    log_fp = out_dir + "split_paired_reads/" + sampleID + "_log.txt"
    log = open(log_fp, 'w')
    subprocess.call(["pefcon", "-1", in_fp, "-o", out_pre, "-z", "--split"],
                    stderr = log)
    out1_fp = out_pre + "_1.fastq"
    out2_fp = out_pre + "_2.fastq"
    gzip1 = subprocess.Popen(["gzip", out1_fp])
    gzip2 = subprocess.Popen(["gzip", out2_fp])
    gzip1.wait()
    gzip2.wait()
    return out1_fp + ".gz", out2_fp + ".gz"


if __name__ == "__main__":
    main()
