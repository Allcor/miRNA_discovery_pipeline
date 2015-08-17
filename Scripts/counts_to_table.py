#!/usr/bin/python

from sys import argv
import operator

def read_fasta(fasta_file):
    with open(fasta_file,'r') as infile:
        fasta_dic = {}
        for line in infile:
            if line[0] == ">":
                name = line.strip().split('_')[0][1:]
                count = line.strip().split('_')[1][1:]
                seq = ""
                fasta_dic[name] = ([count,seq])
            else:
                fasta_dic[name][1] += line.strip()
    return fasta_dic

def make_table(outfile,data):
    with open(outfile,'w') as outfile:
        sorted_data = sorted(data.items(), key=lambda x: int(x[1][0]) , reverse = True)
        for key in sorted_data:
            outfile.write(key[0]+"\t"+data[key[0]][0]+"\t"+data[key[0]][1]+"\n")

def main(arguments):
    fasta_file = arguments[0]
    outfile = arguments[1]
    data = read_fasta(fasta_file)
    make_table(outfile,data)

if __name__ == '__main__':
    main(argv[1:])
