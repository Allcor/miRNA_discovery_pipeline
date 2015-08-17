"""
7-4-2015
    made script, it supposed to translate fastq to fasta, in the format used by mirdeep (with counts)
"""

from sys import (argv,stdout)
from operator import itemgetter
import textwrap

"""READ/WRITE"""
def read_fastq(infile,report):
    counts = {}
    file_type = infile.split('.')[-1]
    if report:
        update_bash_line("loading "+file_type+" file")
    if file_type == "fastq":
        linemarks = ['@','+']
    elif file_type == "fasta":
        linemarks = ['>','>']
    infile = open(infile,'r')
    infile = infile.readlines()
    leng = len(infile)
    is_seq = False
    seq = ""
    if infile[0][0] == linemarks[0]:
        start = infile[0].split(":")[0]
    else:
        print ("can't initialise title convention")
    #going trough file
    for i in range(leng):
        line = infile[i].strip()
        if i % 1000 == 0 and report:
            update_bash_line("logging seqenses: "+str(i)+" / "+str(leng))
        #logging seq
        if line[:len(start)] == start: #indicates title of sequence, start sequence logging
            is_seq = True
            try:
                counts[seq] += 1
            except KeyError:
                counts[seq] = 1
            seq = ""
        elif line[0] == linemarks[1]: #indicates end of sequense
            is_seq = False
        elif is_seq: #logging whole sequence
            seq += line
    #last sequence
    try:
        counts[seq] += 1
    except KeyError:
        counts[seq] = 1
    #for console output:
    if report:
        update_bash_line("logging seqenses: "+str(i)+" / "+str(leng))
    print ("\tDone")
    return counts

def write_fasta(d,report,outfile,name = ""):
    if report:
        update_bash_line("preparing to write file")
    if name == "":
        name = ".".join(outfile.split(".")[:-1])
    outfile = open(outfile,'w')
    count_length = len(str(len(d.keys())))
    reads = sort_on_reads(d)
    counter = 0
    for i in range(len(reads)):
        if i % 100 == 0 and report:
            update_bash_line("writing seqenses: "+str(i)+" / "+str(len(reads)))
        count = reads[i][0]
        seq = reads[i][1]
        counter += 1
        line = ">"+name+str(counter).zfill(count_length)+"_x"+str(count)+"\n"
        outfile.write(line)
        outfile.write("\n".join(textwrap.wrap(seq, width=60))+'\n')
    #for console output:
    if report:
        update_bash_line("writing seqenses: "+str(i)+" / "+str(len(reads)))
    print ("\tDone")

"""TOOLS"""
def sort_on_reads(reads):
    sorted_reads = []
    for key in reads.keys():
        sorted_reads.append((reads[key],key))
    sorted_reads.sort(key=itemgetter(0), reverse=True)
    return sorted_reads

def update_bash_line(string):
    stdout.write("\r"+string)
    stdout.flush()

def help():
    print ("help for mirdeep file maker")
    print ("makes mirdeep format fasta files from fasta and fastq")
    print ("---")
    print ("options:")
    print ("-h or -help")
    print ("-i [file]:\tinput file")
    print ("-o [file]:\toutput file")
    print ("-n [name]:\tsequence names start with this")
    print ("-s :\tsilent mode")

""" MAIN PART"""
def main(arguments):
    if arguments == []:
        help()
    else:
        if arguments[0].split('.')[-1] == "fastq" and arguments[1].split('.')[-1] == "fasta":
            in_file = [arguments[0]]
            out_file = [arguments[1]]
        named = False
        report = True
        for i in range(len(arguments)):
            if arguments[i] == "-i":
                in_file = arguments[i+1].strip().split(',')
            elif arguments[i] == "-o":
                out_file = arguments[i+1].strip().split(',')
            elif arguments[i] == "-h" or arguments[i] == "-help":
                help()
            elif arguments[i] == "-n":
                name = arguments[i+1].strip()
                named = True
            elif arguments[i] == "-s":
                report = False
        for i in range(len(in_file)):
            print ("\nmaking "+in_file[i]+" into "+out_file[i]+"\n")
            #count sequenses
            count_dict = read_fastq(in_file[i],report)
            #write outfile
            if named:
                write_fasta(count_dict,out_file[i],report,name)
            else:
                write_fasta(count_dict,out_file[i],report)

if __name__ == "__main__":
    main(argv[1:])
