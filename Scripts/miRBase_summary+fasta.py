#! /usr/bin/python

'''
This script takes the tables of miRNAs identified with miRBase
and puts them into one long list and also makes a fasta file
out of them that can be used to align/compare with. This
script collapses samples (output has no duplicates) and
numbers of reads are converted into percentages (from the
corresponding sample(s)). Percentages are calculated to mean
percentages if a sequence is present in multiple samples.
Input = miRBase_found_#.csv (where # is all the samples)
Output = miRBase_found-all.csv + miRBase_found-all.fasta
'''

from sys import argv
import os

### argvs to be used are:
# 1: folder with the files (miRBase_found_#.csv)
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
    inputfiles = [read_folder + csv for csv in os.listdir(read_folder) if csv.startswith("miRBase_found_")]
    ## automatically detect files to read from the supplied
    ## folder
    
    writedata = {}
    for inputfile in inputfiles:
    ## loop over each file from the list
        writelist = []
        ## keep a list of the data per file
        total_sequences = 0
        ## keep track of the number of sequences
        with open(inputfile, 'r') as tables:
        ## open each file
            for line in tables:
                line = line.strip('\n').split('\t')
                name = line[0]
                number = int(line[1])
                sequence = line[2]
                writelist.append([name, number, sequence])
                total_sequences += int(number)
                ## read each file and note name, number and 
                ## sequence and append these to the list
                
            print(total_sequences, " in ", inputfile)
                
            for lists in writelist:
                lists[1] = float(lists[1] / total_sequences * 100)
                
                if lists[2] not in writedata:
                    n = 1
                    # keep an 'n' to make a weighted average
                    writedata[lists[2]] = [lists[0], lists[1], n]
                else:
                    n = writedata[lists[2]][2] + 1
                    average = (writedata[lists[2]][1] * n + lists[1]) / (n + 1)
                    writedata[lists[2]] = [writedata[lists[2]][0], average, n]
    
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
        csv.write("Reference\tOccurrence (%)\tSequence\n")
        for keys in data.keys():
            sortlist.append([data[keys][0], data[keys][1], keys])

        sortlist.sort(key=lambda x: x[1], reverse = True)
        ## I want to have the miRNAs sorted by their occurence
        ## and descending (from high to low)
        
        for lists in sortlist:
            csv.write(str(lists[0]) + '\t' + '{0:.5f}'.format(lists[1]) + '\t' + lists[2].replace('T', 'U') + '\n')

    with open(outputfasta, 'w') as fasta:
        for lists in sortlist:
            fasta.write('>' + lists[0] + '_x' + "{0:.5f}%".format(lists[1]) + '\n' + lists[2] + '\n')
        
    print("The list and fasta file are ready!")

def main(argv):
    '''
    This function simply serves to combine all previous
    functions into one working script.
    '''
    write_output(read_lists(argv[0]), argv[0]+argv[1], argv[0]+argv[2])
    
if __name__ == "__main__":
    main(argv[1:])
