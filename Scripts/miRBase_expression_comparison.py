#! /usr/bin/python

'''
Compares expression of miRNAs from different genotypes.
miRBase version

Input = miRBase_found_#.csv (where # is all the samples)
Output = miRBase_expression_comparison.csv
'''


from sys import argv
import os

### argvs to be used are:
# 1: folder with the files (miRBase_found_#.csv)
# 2: output .csv (miRBase_expression_comparison.csv)

def read_files(read_folder, outputfile):
    inputfiles = [read_folder + csv for csv in os.listdir(read_folder) if csv.startswith("miRBase_found_")]

    inputfiles.sort()
    
    carrot_A = [inputfiles[0], inputfiles[3]]
    carrot_B = [inputfiles[1], inputfiles[4]]
    carrot_C = [inputfiles[2], inputfiles[5]]

    total_A = count_total_sequences(carrot_A)
    total_B = count_total_sequences(carrot_B)
    total_C = count_total_sequences(carrot_C)

    print("A:", total_A, "sequences")
    print("B:", total_B, "sequences")
    print("C:", total_C, "sequences")

    total_seqs_dic = {}
            
    seqs_A = collect_sequences(carrot_A, total_A)
    seqs_B = collect_sequences(carrot_B, total_B)
    seqs_C = collect_sequences(carrot_C, total_C)
    
    total_seqs_dic = seqs_A.copy()
    total_seqs_dic.update(seqs_B)
    total_seqs_dic.update(seqs_C)

    header = "Source\tA\tB\tC\tSequence\n"
    
    with open(outputfile, 'w') as outp:
        outp.write(header)
        for sequence in total_seqs_dic.keys():
            if sequence in seqs_A:
                outp.write(seqs_A[sequence][1] + '\t')
                outp.write(str('%.5f' % seqs_A[sequence][0]) + '\t')
            elif sequence in seqs_B:
                outp.write(seqs_B[sequence][1] + '\t-\t')
            elif sequence in seqs_C:
                outp.write(seqs_C[sequence][1] + '\t-\t')
            else:
                outp.write('-\t')
            if sequence in seqs_B:
                outp.write(str('%.5f' % seqs_B[sequence][0]) + '\t')
            else:
                outp.write('-\t')
            if sequence in seqs_C:
                outp.write(str('%.5f' % seqs_C[sequence][0]) + '\t')
            else:
                outp.write('-\t')
            outp.write(sequence)
                
    print("Unique sequences:", len(total_seqs_dic), '\n')
    
def collect_sequences(inputfiles, total):
    '''
    Collects all the sequences and expression levels for the
    given carrot, from the sequences that were foundn with
    miRBase.
    Returns a dictionary of {sequence: expression, source 
    name} for each sequence present.
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
                sequence = line[2].replace('T', 'U')
                source_name = line[0]
                expression = int(line[1]) / total * 100
                if sequence in seq_dictio:
                    seq_dictio[sequence][0] = (seq_dictio[sequence][0] + expression) / 2
                else:
                    seq_dictio[sequence] = [expression, source_name]
    return seq_dictio
    
def count_total_sequences(inputfiles):
    '''
    Count the total number of sequences for each carrot.
    Returns the number of sequences (int).
    inputfiles = a list of files from which to read/count/
     sum all the sequences.
    '''
    total_number = 0
    for inputfile in inputfiles:
        with open(inputfile, 'r') as tables:
            tables.readline()
            ## skip the first line (header)
            for line in tables:
                line = line.strip('\n').split('\t')
                number = int(line[1])
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
