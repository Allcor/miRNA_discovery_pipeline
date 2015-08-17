#! /usr/bin/python

from sys import argv

script, inputfile, outputfile = argv

def miRDP_to_fasta_mature(inputfile, outputfile):
    with open(inputfile, 'r') as readfile:
        with open(outputfile, 'w') as writefile:
            for line in readfile.readlines():
                name = line.split()[2]
                sequence = line.split()[6]
                writefile.write('>' + name + '\n' + sequence + '\n')
    
miRDP_to_fasta_mature(inputfile, outputfile)
