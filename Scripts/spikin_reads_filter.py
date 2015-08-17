#!/usr/bin/python
'''
we have a bunch of RNA-seq data (fastq file) and want to filter out known sequences (fasta).

4-3-2015
added a couple of lines to prevent filtering of negetive strand mapped reads, this can be usefull when the database used only has relevent invormation on the + strand.

80% rule thing, lenth of mapped reads need to be 80% of the spike in or it is not counted.

23-4-2015
made it faster by not using .keys() ...
also added a size cutoff
30-4-2015
options properly inplemented
'''

SPIKES_NSPK = ['NSPK_RNA_1', 'NSPK_RNA_2', 'NSPK_RNA_3', 'NSPK_RNA_4', 'NSPK_RNA_5', 'NSPK_RNA_6', 'NSPK_RNA_7', 'NSPK_RNA_8', 'NSPK_RNA_9', 'NSPK_RNA_10', 'NSPK_RNA_11', 'NSPK_RNA_12', 'NSPK_RNA_13', 'NSPK_RNA_14', 'NSPK_RNA_15', 'NSPK_RNA_16', 'NSPK_RNA_17', 'NSPK_RNA_18','NSPK_RNA_19','NSPK_RNA_20']
SPIKES_SSPK = ['SSPK_RNA_10', 'SSPK_RNA_16', 'SSPK_RNA_19', 'SSPK_RNA_22', 'SSPK_RNA_25', 'SSPK_RNA_28', 'SSPK_RNA_34', 'SSPK_RNA_40', 'SSPK_RNA_50', 'SSPK_RNA_60', 'SSPK_RNA_70']

from sys import argv,stdout
import os

def log_sam(sam, silent):
    sequences = {}
    minimal_length = {}
    counter = 0
    for line in sam:
        counter += 1
        if counter % 100 == 0 and silent != True:
            update_bash_line(str(counter)+" lines read,")
        splitline = line.strip().split("\t")
            #no negetive strand mapped reads allowed
        if splitline[1] != "4" and line[0] != "@" and splitline[1] != "16":
            if len(splitline[9]) >= minimal_length[splitline[2]]:
                sequences[splitline[9]] = line
        elif line[:3] == "@SQ":
            size = int(line.split(":")[-1].strip())
            if size <= 40:
                min_size = round(0.8 * size)
            else:
                min_size = 32
            minimal_length[line.split(":")[1].split("\t")[0]] = min_size
    if silent != True:
        stdout.write(" done\n")
    return sequences

def log_fasta(fasta):
    sequences = {}
    seq = ""
    for line in fasta:
        if line[0] != ">":
            seq += line.strip()
        else:
            if seq != "":
                sequences[seq] = line
            seq = ""
    return sequences

def update_bash_line(string):
    stdout.write("\r"+string)
    stdout.flush()

def filter_fastq(fastq, targets, outfile, silent, cutoff = [8,'all']):
    end = False
    filtered = {}
    keys = []
    label = ""
    seq = ""
    quality = ""
    lower_count = 0
    upper_count = 0
    counter = 0
    with open(outfile,"a+") as append_line:
        while not end:
            counter += 1
            if counter % 100 == 0 and silent != True:
                update_bash_line(str(counter)+" fastq sequences checked,")
            line = fastq.readline()
            if line == "":
                end = True
            elif line[0] == "@":
                seq = fastq.readline().strip()
                seq_len = len(seq)
                label = fastq.readline()[1:].strip()
                quality = fastq.readline().strip()
                write = True
                if seq_len < int(cutoff[0]):
                    write = False
                    lower_count += 1
                elif cutoff[1] != 'all':
                    if seq_len > int(cutoff[1]):
                        write = False
                        upper_count += 1
                if seq in targets:
                    hit = False
                    spike_name = targets[seq].strip().split("\t")[2]
                    if spike_name in filtered:
                        filtered[spike_name] += 1
                    else:
                        filtered[spike_name] = 1
                if write:
                    append_line.write("@"+label+"\n")
                    append_line.write(seq+"\n")
                    append_line.write("+"+label+"\n")
                    append_line.write(quality+"\n")
            else:
                print (line)
    if silent != True:
        stdout.write("\tdone\n")
    return filtered,lower_count,upper_count
    
def report(filtered,report_name,in_name):
    directory = '/'.join(report_name.split('/')[:-1])
    if not os.path.exists(directory):
        os.makedirs(directory)
    if os.path.isfile(report_name) == False:
        with open(report_name,"a+") as report_file:
            report_file.write('\t'+'\t'.join(SPIKES_NSPK)+"\t|\t"+'\t'.join(SPIKES_SSPK)+"\n")
    for colomn in SPIKES_NSPK:
        if colomn not in filtered:
            filtered[colomn] = 0
    for colomn in SPIKES_SSPK:
        if colomn not in filtered:
            filtered[colomn] = 0
    with open(report_name,"a+") as report_file:
        report_file.write(in_name+'\t'+'\t'.join([str(filtered[colomn]) for colomn in SPIKES_NSPK])+"\t|\t")
        report_file.write('\t'.join([str(filtered[colomn]) for colomn in SPIKES_SSPK])+'\n')

def removed_stuff_log(file_name,l_count,u_count,in_name):
    if os.path.isfile(file_name) == False:
        with open(file_name,'w') as document:
            document.write("file_name\treads_removed_lower\treads_removed_upper\n")
    with open(file_name,'a+') as document:
        document.write(in_name+'\t'+str(l_count)+'\t'+str(u_count)+'\n')
 
def main(arguments):
    report_name = False
    length_cutoff = False
    filter_doc = False
    silent = False
    #### options
    for i in range(len(arguments)):
        argument = arguments[i]
        if argument == "-l":
            try:
                length_cutoff = arguments[i+1].split(',')
            except KeyError:
                None
        if argument == "-i":
            infile = arguments[i+1]
        if argument == "-q":
            fastq = arguments[i+1]
            fastq_file = open(fastq, 'r')
        if argument == "-r":
            report_name = arguments[i+1]
        if argument == "-o":
            new_fastq = arguments[i+1]
        if argument == "-p":
            filter_doc = arguments[i+1]
        if argument == "-s":
            silent = True
    ######
    if infile[-4:] == ".sam":
        sam_file = open(infile, 'r')
        seq_to_filter = log_sam(sam_file,silent)
    elif infile[-6:] == '.fasta':
        fasta_file = open(infile, 'r')
        seq_to_filter = log_seqs(fasta_file.readlines())
    ######
    print ("filtering "+fastq.split('/')[-1])
    file_number = infile.split('/')[-1].split('_')[1]
    if length_cutoff != False:
        filtered,size_lower_count,size_upper_count = filter_fastq(fastq_file, seq_to_filter, new_fastq , silent, length_cutoff)
    else:
        filtered,size_lower_count,size_upper_count = filter_fastq(fastq_file, seq_to_filter, new_fastq, silent)
    if report_name != False:
        report(filtered,report_name,file_number)
    all_filtered = 0
    for i in filtered:
        all_filtered += int(filtered[i])
    if filter_doc != False:
        removed_stuff_log(filter_doc,size_lower_count,size_upper_count,file_number)
    print ("")
    print ("sequences filtered out:")
    print ("total spikes: ", all_filtered)
    print ("lower size restriction: ", size_lower_count)
    print ("upper size restriction: ", size_upper_count)

if __name__ == "__main__":
    main(argv[1:])
