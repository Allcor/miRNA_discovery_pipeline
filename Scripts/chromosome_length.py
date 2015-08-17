#!/usr/bin/python

"""
gous trough the genome fasta and extracts the length.
"""


from sys import argv

script,genome,output = argv

with open(genome,'r') as infile:
    with open(output,'w') as outfile:
        seq_info = ["",[]]
        for line in infile:
            if line[0] == ">":
                if seq_info[0] != "":
                    name = seq_info[0][1:].split()[0]
                    length = len("".join(seq_info[1]))
                    outfile.write(name+"\t"+str(length)+"\n")
                seq_info = [line,[]]
            else:
                seq_info[1].append(line.strip())

