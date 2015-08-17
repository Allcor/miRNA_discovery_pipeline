#! /usr/bin/python

'''
This script takes the tables of miRNAs identified with 
miRDeep-P and puts them into one long list and also makes a 
fasta file out of them that can be used to align/compare 
with. This script collapses samples (output has no
duplicates) and numbers of reads are converted into
percentages (from the corresponding sample(s)). Percentages 
are calculated to mean percentages if a sequence is present 
in multiple samples.
Input = Sample_#_mature-miRBase.csv
        Sample_#_mature-notinmiRBase.csv
         (where # is all the samples)
Output = miRDP_found-all.csv + miRDP_found-all.fasta
'''

from sys import argv
import os

### argvs to be used are:
# 1: folder with the files (Sample_#_mature-*.csv)
# 2: output .csv
# 3: output .fasta

def read_lists(read_folder):
    '''
    This function reads the separate files with miRBase hits.
    It returns a dictionary of data to be written by the next
    function.
    Input = miRBase_found_#.csv (files)
    Output = dictionary of data to be used
    '''
    inputfiles_miRBase = [read_folder + csv for csv in os.listdir(read_folder) if csv.endswith("-miRBase.csv")]
    inputfiles_notIn = [read_folder + csv for csv in os.listdir(read_folder) if csv.endswith("_notinmiRBase.csv")]
    ## automatically detect files to read from the supplied
    ## folder
    
    writedata = {}
    total_sequences = 0
    ## keep track of the number of sequences
    for inputfile in inputfiles_miRBase:
        with open(inputfile, 'r') as tables:
            tables.readline()
            ## skip the first line (header)
            for line in tables:
                line = line.strip('\n').split('\t')
                number = int(line[10])
                total_sequences += number
            print(total_sequences, " in ", inputfile)
                
    for inputfile in inputfiles_notIn:
        with open(inputfile, 'r') as tables:
            tables.readline()
            ## skip the first line (header)
            for line in tables:
                line = line.strip('\n').split('\t')
                number = int(line[8])
                total_sequences += number
            print(total_sequences, " in ", inputfile)

    for inputfile in inputfiles_miRBase:
    ## loop over each file from the list
        writelist = []
        ## keep a list of the data per file
        with open(inputfile, 'r') as tables:
        ## open each file
            tables.readline()
            ## skip the first line (header)
            for line in tables:
                line = line.strip('\n').split('\t')
                source_name = line[0]
                name = line[1]
                prec_len = line[2]
                MFEI = line[5]
                sequence = line[6]
                mismatches = line[7]
                length = line[8]
                mat_seq_arm = line[9]
                number = int(line[10])
                writelist.append([name, number, sequence])
                ## read each file and note name, number and 
                ## sequence and append these to the list
                
            for lists in writelist:
                lists[1] = float(lists[1] / total_sequences * 100)
                
                if lists[2] not in writedata:
                    writedata[lists[2]] = [lists[0], lists[1]]
                else:
                    sum = (writedata[lists[2]][1] + lists[1])
                    writedata[lists[2]] = [writedata[lists[2]][0], sum]

    for inputfile in inputfiles_notIn:
    ## loop over each file from the list
        writelist = []
        ## keep a list of the data per file
        with open(inputfile, 'r') as tables:
        ## open each file
            tables.readline()
            ## skip the first line (header)
            for line in tables:
                line = line.strip('\n').split('\t')
                name = line[0]
                prec_len = line[1]
                MFEI = line[4]
                sequence = line[5]
                length = line[6]
                mat_seq_arm = line[7]
                number = int(line[8])
                writelist.append([name, number, sequence])
                ## read each file and note name, number and 
                ## sequence and append these to the list
                
            for lists in writelist:
                lists[1] = float(lists[1] / total_sequences * 100)
                
#                if lists[2] not in writedata:
#                    n = 1
#                    # keep an 'n' to make a weighted average
#                    writedata[lists[2]] = [lists[0], lists[1], n]
#                else:
#                    n = writedata[lists[2]][2] + 1
#                    average = (writedata[lists[2]][1] * n + lists[1]) / (n + 1)
#                    writedata[lists[2]] = [writedata[lists[2]][0], average, n]

                if lists[2] not in writedata:
                    writedata[lists[2]] = [lists[0], lists[1]]
                else:
                    sum = (writedata[lists[2]][1] + lists[1])
                    writedata[lists[2]] = [writedata[lists[2]][0], sum]

    print("Total number of sequences identified as miRNA:", total_sequences)
    return writedata
    
def write_output(data, outputcsv, outputfasta):
    '''
    This function writes all the required data to a new,
    collapsed .csv file and to a fasta file. It takes the
    dictionary created by the function above.
    Input = dictionary of data to be used
    Output = miRBase_found_all.csv (file) +
     miRBase_found_all.fasta (file)
    '''
    sortlist = []
    with open(outputcsv, 'w') as csv:
        csv.write("Name\tOccurrence (%)\tSequence\n")
        for keys in data.keys():
            sortlist.append([data[keys][0], data[keys][1], keys])

        sortlist.sort(key=lambda x: x[1], reverse = True)
        ## I want to have the miRNAs sorted by their occurence
        ## and descending (from high to low)
        
        for lists in sortlist:
            csv.write(str(lists[0]) + '\t' + '{0:.5f}'.format(lists[1]) + '\t' + lists[2] + '\n')

    with open(outputfasta, 'w') as fasta:
        for lists in sortlist:
            fasta.write('>' + lists[0] + '_x' + "{0:.5f}%".format(lists[1]) + '\n' + lists[2].replace('U', 'T').replace('u', 't') + '\n')
        
    print("The list and fasta file are ready!")

def main(argv):
    '''
    This function simply serves to combine all previous
    functions into one working script.
    '''
    write_output(read_lists(argv[0]), argv[0]+argv[1], argv[0]+argv[2])
    
if __name__ == "__main__":
    main(argv[1:])
