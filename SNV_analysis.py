#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 17:42:06 2022

@author: seob
"""

import pandas as pd
import pysam
import pickle
import os


def get_reference_position(bam_file_list, save_dir):
    '''
    Input : bam_file_list: list of the bamfile path
    output: save the pickle data in the save_dir
        pickle data is the list of dict of each read data in bam file
        dict['chr'] : chormosome of the read
        dict['start'] : start location of the reference sequence
        dict['end'] : end location of the reference sequence
        dict['seq'] : sequence of the aligned sequence
    '''
    name_list = []
    for index, file in enumerate(bam_file_list):
        file_name = file.split('/')[-1].split('.')[0]
        name_list.append(file_name)
        samfile = pysam.AlignmentFile(file)
        read_list = []
        for read in samfile.fetch():
            if read.is_forward:
                read_dict = {}
                end = read.reference_end
                data = read.get_aligned_pairs()
                seq = read.get_forward_sequence()
                min_pair = -1
                for paired in data:
                    if paired[1] != None and paired[0] != None:
                        if min_pair == -1:
                            min_pair = paired[1]
                            min_start = paired[0]
                        max_pair = paired[1]
                        end = paired[0]
                read_dict['chr'] = read.reference_id + 1
                read_dict['start'] = min_pair
                read_dict['end'] = max_pair
                read_dict['seq'] = seq[min_start:end]
                read_list.append(read_dict)

        with open(os.path.join(save_dir, name_list[index] + '.pickle'),'wb') as f:
            pickle.dump(read_list, f)


bam_file_list = ["gDNA.uniqMap.bam"]
pickle_save_dir = "pickle"
get_reference_position(bam_file_list, pickle_save_dir)
target_seq_data = "target_seqs.csv"
target_seq_df = pd.read_csv(target_seq_data)

data_list = os.listdir(pickle_save_dir)
name_list = [data.split('.')[0] for data in data_list]
data_list = [os.path.join(pickle_save_dir, data) for data in data_list]
final_list = []

for name in name_list:
    print(name)
    with open(os.path.join(pickle_save_dir, name + '.pickle'), "rb") as fr:
        data = pickle.load(fr)
    aligned_df = pd.DataFrame(data)
    snp_list = []
    snv_list = []
    gene_list = []
    ref_gene_dict = {}
    for index, row in target_seq_df.iterrows():
        chrom = row['Chromosome']
        if chrom == 'X':
            chrom = 23
        if chrom != 'M':
            chrom = int(chrom)
        chrom_df = aligned_df[aligned_df['chr'] == chrom]
        position = row['Start_position']
        alt = row['Tumor_Seq_Allele1'] # Variant allele sequence
        ref = row['Reference_Allele'] # Reference sequence
        gene = row['Hugo_Symbol'] # Gene name

        if gene not in ref_gene_dict.keys():
            ref_gene_dict[gene] = [position]
            if ref != '-' and alt != '-':
                if len(ref) == 1 and len(alt) == 1:
                    pos_df = chrom_df[chrom_df['start'] < position]
                    pos_df = pos_df[pos_df['end'] > position]
                    if len(pos_df) != 0:
                        pcr_start_list = []
                        total_read_count = 0
                        snv_read_count = 0
                        ref_read_count = 0
                        for data_idx, data_row in pos_df.iterrows():
                            seq_data = data_row['seq']
                            pcr_start_list.append(data_row['start'])
                            total_read_count += 1
                            location = position - data_row['start']
                            if len(seq_data) - 1 > location:
                                sequence = seq_data[location - 1]
                                if sequence == alt:
                                    snv_list.append(gene + '_' + str(position) + '_' + ref + '_' + alt)
                                    gene_list.append(gene)
                                    snv_read_count += 1
                                if sequence == ref:
                                    snp_list.append({'gene': gene, 'pos': position})
                                    ref_read_count += 1
                        pcr_read_count = len(pcr_start_list) - len(set(pcr_start_list))
                        final_list.append(
                            {'cell': name, 'gene': gene, 'chr': chrom, 'pos': position, 'seq': sequence, 'ref': ref,
                             'alt': alt, 'total_count': total_read_count,
                             'PCR_read_count': pcr_read_count, 'SNV_read_count': snv_read_count,
                             'ref_read_count': ref_read_count})

        elif position not in ref_gene_dict[gene]:
            ref_gene_dict[gene].append(position)
            if ref != '-' and alt != '-':
                if len(ref) == 1 and len(alt) == 1:
                    pos_df = chrom_df[chrom_df['start'] < position]
                    pos_df = pos_df[pos_df['end'] > position]
                    if len(pos_df) != 0:
                        pcr_start_list = []
                        total_read_count = 0
                        snv_read_count = 0
                        ref_read_count = 0
                        for data_idx, data_row in pos_df.iterrows():
                            seq_data = data_row['seq']
                            pcr_start_list.append(data_row['start'])
                            total_read_count += 1
                            location = position - data_row['start']
                            if len(seq_data) - 1 > location:
                                sequence = seq_data[location - 1]
                                if sequence == alt:
                                    snv_list.append(gene + '_' + str(position) + '_' + ref + '_' + alt)
                                    gene_list.append(gene)
                                    snv_read_count += 1
                                if sequence == ref:
                                    snp_list.append({'gene': gene, 'pos': position})
                                    ref_read_count += 1
                        pcr_read_count = len(pcr_start_list) - len(set(pcr_start_list))
                        final_list.append(
                            {'cell': name, 'gene': gene, 'chr': chrom, 'pos': position, 'seq': sequence, 'ref': ref,
                             'alt': alt, 'total_count': total_read_count,
                             'PCR_read_count': pcr_read_count, 'SNV_read_count': snv_read_count,
                             'ref_read_count': ref_read_count})
    print(len(set(snv_list)))
pandas_data = pd.DataFrame(final_list)
result_dir = 'result'
pandas_data.to_csv(os.path.join(result_dir, 'SNV_data.csv'), index = None)