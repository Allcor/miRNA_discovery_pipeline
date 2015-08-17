#! /usr/bin/python

'''
Adjusted version of 'miRDP_summary+fasta.py' to make one
long list with all (relevant) details.
'''

from sys import argv
import os

### argvs to be used are:
# 1: folder with the files (Sample_#_mature-*.csv)
# 2: output .csv

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
                mfe = line[3]
                MFEI = line[5]
                sequence = line[6]
                mismatches = line[7]
                length = line[8]
                mat_seq_arm = line[9]
                number = int(line[10])
                writelist.append([source_name, name, prec_len, mfe, MFEI, sequence, length, mat_seq_arm])
                ## read each file and note name, number and 
                ## sequence and append these to the list
                
            for lists in writelist:
                if lists[5] not in writedata:
                    writedata[lists[5]] = [lists[0], lists[1], lists[2], lists[3], lists[4], lists[6], lists[7]]
                else:
                    pass

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
                mfe = line[3]
                MFEI = line[4]
                sequence = line[5]
                length = line[6]
                mat_seq_arm = line[7]
                number = int(line[8])
                writelist.append(['-', name, prec_len, mfe, MFEI, sequence, length, mat_seq_arm])
                ## read each file and note name, number and 
                ## sequence and append these to the list
                
            for lists in writelist:
                if lists[5] not in writedata:
                    writedata[lists[5]] = [lists[0], lists[1], lists[2], lists[3], lists[4], lists[6], lists[7]]
                else:
                    pass

    print("Total number of sequences identified as miRNA:", total_sequences)
    return writedata
    
def write_output(data, outputcsv):
    '''
    This function writes all the required data to a new,
    collapsed .csv file. It takes the
    dictionary created by the function above.
    Input = dictionary of data to be used
    Output = miRBase_found_all.csv (file) +
    '''
    sortlist = []
    with open(outputcsv, 'w') as csv:
        csv.write("Source name\tName\tPrec Len\tMin Free Energy\tMFEI\tSequence\tLength\tMat Seq Arm\n")
        for keys in data.keys():
            sortlist.append([data[keys][0], data[keys][1], data[keys][2], data[keys][3], data[keys][4], data[keys][5], keys, data[keys][6]])

        sortlist.sort(key=lambda x: x[0], reverse = True)
        
        for lists in sortlist:
            csv.write('\t'.join(lists) + '\n')
        
    print("The list is ready!")

def main(argv):
    '''
    This function simply serves to combine all previous
    functions into one working script.
    '''
    write_output(read_lists(argv[0]), argv[0]+argv[1])
    
if __name__ == "__main__":
    main(argv[1:])
