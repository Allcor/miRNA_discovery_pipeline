#! /usr/bin/python

'''
Compares expression of miRNAs from different genotypes.
miRDP version

Input = Sample_#_mature-miRBase.csv
        Sample_#_mature-notinmiRBase.csv
         (where # is all the samples)
Output = miRDP_expression_comparison.csv
'''

from sys import argv
import os

### argvs to be used are:
# 1: folder with the files (Sample_#_mature-*.csv)
# 2: output .csv (miRDP_expression_comparison.csv)

def read_files(read_folder, outputfile):
    inputfiles_miRBase = [read_folder + csv for csv in os.listdir(read_folder) if csv.endswith("-miRBase.csv")]
    inputfiles_notIn = [read_folder + csv for csv in os.listdir(read_folder) if csv.endswith("_notinmiRBase.csv")]
    
    inputfiles_miRBase.sort()
    inputfiles_notIn.sort()
    
    carrot_A_miRBase = [inputfiles_miRBase[0], inputfiles_miRBase[3]]
    carrot_B_miRBase = [inputfiles_miRBase[1], inputfiles_miRBase[4]]
    carrot_C_miRBase = [inputfiles_miRBase[2], inputfiles_miRBase[5]]
    
    carrot_A_notIn = [inputfiles_notIn[0], inputfiles_notIn[3]]
    carrot_B_notIn = [inputfiles_notIn[1], inputfiles_notIn[4]]
    carrot_C_notIn = [inputfiles_notIn[2], inputfiles_notIn[5]]

    total_A = count_total_sequences(carrot_A_miRBase) + count_total_sequences(carrot_A_notIn)
    total_B = count_total_sequences(carrot_B_miRBase) + count_total_sequences(carrot_B_notIn)
    total_C = count_total_sequences(carrot_C_miRBase) + count_total_sequences(carrot_C_notIn)
            
    print("A:", total_A, "sequences")
    print("B:", total_B, "sequences")
    print("C:", total_C, "sequences")

    total_seqs_dic = {}
            
    seqs_A_miRBase = collect_sequences_miRBase(carrot_A_miRBase, total_A)
    seqs_A_notIn = collect_sequences_notIn(carrot_A_notIn, total_A)
    all_seqs_A = seqs_A_miRBase.copy()
    all_seqs_A.update(seqs_A_notIn)

    seqs_B_miRBase = collect_sequences_miRBase(carrot_B_miRBase, total_B)
    seqs_B_notIn = collect_sequences_notIn(carrot_B_notIn, total_B)
    all_seqs_B = seqs_B_miRBase.copy()
    all_seqs_B.update(seqs_B_notIn)
    
    seqs_C_miRBase = collect_sequences_miRBase(carrot_C_miRBase, total_C)
    seqs_C_notIn = collect_sequences_notIn(carrot_C_notIn, total_C)
    all_seqs_C = seqs_C_miRBase.copy()
    all_seqs_C.update(seqs_C_notIn)
    
    total_seqs_dic = all_seqs_A.copy()
    total_seqs_dic.update(all_seqs_B)
    total_seqs_dic.update(all_seqs_C)

    header = "Sequence\tA\tB\tC\n"
    
    with open(outputfile, 'w') as outp:
        outp.write(header)
        for sequence in total_seqs_dic.keys():
            outp.write(sequence + '\t')
            if sequence in all_seqs_A:
                outp.write(str(all_seqs_A[sequence]) + '\t')
            else:
                outp.write('-\t')
            if sequence in all_seqs_B:
                outp.write(str(all_seqs_B[sequence]) + '\t')
            else:
                outp.write('-\t')
            if sequence in all_seqs_C:
                outp.write(str(all_seqs_C[sequence]) + '\n')
            else:
                outp.write('-\n')
                
    print("Unique sequences:", len(total_seqs_dic), '\n')
    
def collect_sequences_miRBase(inputfiles, total):
    '''
    Collects all the sequences and expression levels for the
    given carrot, from the sequences that are in miRBase.
    Returns a dictionary of {sequence: expression} lists for
    each sequence present.
    inputfiles = list of files
    total = number of sequences from a carrot line (int)
    '''
    seq_dictio = {}
    for files in inputfiles:
        with open(files, 'r') as readfile:
            readfile.readline()
            # skip the first line (header)
            for line in readfile:
                line = line.split('\t')
                sequence = line[6]
                expression = int(line[10]) / total * 100
                if sequence in seq_dictio:
                    seq_dictio[sequence] = (seq_dictio[sequence] + expression) / 2
                else:
                    seq_dictio[sequence] = expression
    return seq_dictio
    
def collect_sequences_notIn(inputfiles, total):
    '''
    Collects all the sequences and expression levels for the
    given carrot, from the sequences that are *not* in miRBase.
    Returns a dictionary of {sequence: expression} lists for
    each sequence present.
    inputfiles = list of files
    total = number of sequences from a carrot line (int)
    '''
    seq_dictio = {}
    for files in inputfiles:
        with open(files, 'r') as readfile:
            readfile.readline()
            # skip the first line (header)
            for line in readfile:
                line = line.split('\t')
                sequence = line[5]
                expression = int(line[8]) / total * 100
                if sequence in seq_dictio:
                    seq_dictio[sequence] = (seq_dictio[sequence] + expression) / 2
                else:
                    seq_dictio[sequence] = expression
    return seq_dictio
    
def count_total_sequences(inputfiles):
    '''
    Count the total number of sequences for each carrot.
    Returns the number of sequences (int).
    inputfiles = a list of files from which to read/count/
     sum all the sequences. Includes both the -miRBase.csv 
     files and the _notinmiRBase.csv files.
    '''
    total_number = 0
    for inputfile in inputfiles:
        if inputfile.endswith("-miRBase.csv"):
            with open(inputfile, 'r') as tables:
                tables.readline()
                ## skip the first line (header)
                for line in tables:
                    line = line.strip('\n').split('\t')
                    number = int(line[10])
                    total_number += number
                
        else:
            with open(inputfile, 'r') as tables:
                tables.readline()
                ## skip the first line (header)
                for line in tables:
                    line = line.strip('\n').split('\t')
                    number = int(line[8])
                    total_number += number
                
    return total_number

def main(argv):
    '''
    This function simply serves to combine all previous
    functions into one working script.
    '''
    read_files(argv[0], (argv[0]+argv[1]))

if __name__ == "__main__":
    main(argv[1:])
