#! /usr/bin/python

'''
Analysis of Barozai's miRNAs in our miRNAs identified with
miRBase. Counts hits by bowtie2 and summarises results per 
carrot line and gives a total overview
 (Do all samples contain the 17 miRNAs?).

input = .sam files by alignments of each sample (miRNAs
 identified with miRBase)
output = summary of Barozai's miRNAs found in our miRNAs,
 e.g. "sample 1: xxx sequences
       sample 2: xxx sequences
       sample 3: xxx sequences
       sample 4: xxx sequences
       sample 5: xxx sequences
       sample 6: xxx sequences
       
       A: 17/17
       B: 17/17
       C: 17/17
       TOTAL: All miRNAs were found in all samples."
'''

import os

FASTQ_FOLDER = "Databases/raw/"
SAMPLES =  [s[:-6] for s in os.listdir(FASTQ_FOLDER) if s.endswith(".fastq")]
SAM_FILES = ["./Mapping_results/sample_", "_miRBase-pakibase.sam"]
OUTPUT = "Report/miRNAs-to-pakibase_with-miRBase.txt"

miRNAs_in_sample = {}
sequences_per_sample = {}
miRNAs_in_Base = {}

for sample in SAMPLES:
    names_counted = []
    miRNAs_in_Base[sample] = []
    miRNAs_in_sample[sample] = 0
    sequences_per_sample[sample] = 0
    xsample = open(SAM_FILES[0] + sample + SAM_FILES[1], 'r')
    for line in xsample:
        line = line.strip('\n')
        if line[0] == '@':
            if line[1:3] == 'SQ':
                miRNAs_in_Base[sample].append(line[7:-6])
            else:
                pass
        else:
            sequences_per_sample[sample] += 1
            miRNA_name = line.split('\t')[2]
            if miRNA_name in names_counted:
                pass
            else:
                names_counted.append(miRNA_name)
    miRNAs_in_sample[sample] = names_counted

miRNAs_in_A = list(set(miRNAs_in_sample['1'] + miRNAs_in_sample['4']))
no_miRNAs_in_A = len(miRNAs_in_A)
miRNAs_in_B = list(set(miRNAs_in_sample['2'] + miRNAs_in_sample['5']))
no_miRNAs_in_B = len(miRNAs_in_B)
miRNAs_in_C = list(set(miRNAs_in_sample['3'] + miRNAs_in_sample['6']))
no_miRNAs_in_C = len(miRNAs_in_C)

all_seqs = len(miRNAs_in_Base['1'])

if no_miRNAs_in_A ==all_seqs and no_miRNAs_in_B == all_seqs and no_miRNAs_in_C == all_seqs:
    summary = "All miRNAs found in all samples!"
else:
    summary = """\nA: %i/%i miRNAs found.
B: %i/%i miRNAs found.
C: %i/%i miRNAs found.""" % (no_miRNAs_in_A, all_seqs,
 no_miRNAs_in_B, all_seqs, no_miRNAs_in_C, all_seqs)
    miRNAs_not_found = {}
    for sample in SAMPLES:
        miRNAs_not_found[sample] = []
        for miRNAs in miRNAs_in_Base[sample]:
            if miRNAs not in miRNAs_in_sample[sample]:
                miRNAs_not_found[sample].append(miRNAs)
    
    miRNAs_not_in_A = []
    miRNAs_not_in_B = []
    miRNAs_not_in_C = []
    
    for miRNAs in miRNAs_in_Base['1']:
        if miRNAs not in miRNAs_in_A:
            miRNAs_not_in_A.append(miRNAs)
        if miRNAs not in miRNAs_in_B:
            miRNAs_not_in_B.append(miRNAs)
        if miRNAs not in miRNAs_in_C:
            miRNAs_not_in_C.append(miRNAs)
            
    summary += """\n\nmiRNAs NOT in A: %s
miRNAs NOT in B: %s
miRNAs NOT in C: %s""" % (', '.join(miRNAs_not_in_A), ', '.join(miRNAs_not_in_B), ', '.join(miRNAs_not_in_C))
                
with open(OUTPUT, 'w') as outputfile:
    outputfile.write("Sample 1:\t %i sequences\n" % sequences_per_sample['1'])
    outputfile.write("Sample 2:\t %i sequences\n" % sequences_per_sample['2'])
    outputfile.write("Sample 3:\t %i sequences\n" % sequences_per_sample['3'])
    outputfile.write("Sample 4:\t %i sequences\n" % sequences_per_sample['4'])
    outputfile.write("Sample 5:\t %i sequences\n" % sequences_per_sample['5'])
    outputfile.write("Sample 6:\t %i sequences\n" % sequences_per_sample['6'])
    outputfile.write(summary)
