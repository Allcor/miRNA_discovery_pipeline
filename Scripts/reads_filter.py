#!/usr/bin/python

'''
With the DataQC we made a set of fastq files of our sequencing data, this script is supposed to make a new fastq with mapped (bowtie2 sam file) filtered out. 
- want to start with labeling the reads. reference and MAPQ (MAPing Qualety), it is on position 4(5th entry).
30-3-2015 : report of all labels added to the fastq
31-3-2015 : can also filter (-f)
29-4-2015 : made it more modular, now can map multiple sam files.
'''

import random
from sys import (argv,stdout)

def read_sam(sam_file,progress):
    dic = {}
    sam_name = sam_file.split('/')[-1]
    sam = open(sam_file, "r").readlines()
    leng = len(sam)
    for i in range(leng):
        line = sam[i]
        if i%100 == 0 and progress:
            update_bash_line("reading " + sam_name + " file: "+str(i)+" / "+str(leng))
        if line[0] != "@":
            line = line.strip().split()
            if line[1] != "4":
                dic[line[0].split("|")[0].strip()] = (line[2],line[4])
    if progress:
        update_bash_line("reading " + sam_name + " file: "+str(i+1)+" / "+str(leng))
        print ("\tDone")
    return dic

def stats_and_write(fastq,output,progress,filters = []):
    leng = len(fastq)
    out = open(output,"w")
    count = 0
    new_count = 0
    flag_lib = {}
    skip = True
    start = True
    for i in range(leng):
        line = fastq[i].strip()
        if start == True and line[0] == '@':
            start = fastq[0].split(":")[0]
        else: ("skipping : "+line)
        if line[:len(start)] == start:
            #count total sequenses
            count += 1
            if i % 100 == 0 and progress:
                update_bash_line("writing new fastq file: "+str(i+1)+" / "+str(leng))
            #count 'best' label
            labels = line.split("|")[1:]
            skip = False
            if labels != []:
                if len(labels) > 1: 
                    target_score = 255
                    multiple = 1
                    for label in labels:
                        score = int(label.split("-")[-1].strip())
                        if score < target_score:
                            label_to_log = label
                            target_score = score
                        elif score == target_score: # if there is no clear winner, one is picked at random
                            if type(label_to_log) == "list":
                                label_to_log += [label]
                            elif type(label_to_log) == "str":
                                label_to_log = [label_to_log,label]
                    if type(label_to_log) == "list":
                        label_to_log = random.choice(label_to_log)
                else:
                    multiple = 0
                    label_to_log = labels[0]
                if len(label_to_log.split("-")) != 2: #if '-' is in the name part of the label, replaced by '_'
                    label_to_log = "_".join(label_to_log.strip().split("-")[:-1])+" - "+label_to_log.split("-")[-1]
                # check if needs to be skipped
                if label_to_log.split("-")[0].strip() in filters:
                    skip = True
                    flag = label_to_log.split("-")[0].strip()
                    try:
                        flag_lib[flag+", deleted"][multiple] += 1
                    except KeyError:
                        flag_lib[flag+", deleted"] = [0,0]
                        flag_lib[flag+", deleted"][multiple] += 1
                if not skip:
                    # add to lib
                    flag = label_to_log.split("-")[0].strip()
                    try:
                        flag_lib[flag][multiple] += 1
                    except KeyError:
                        flag_lib[flag] = [0,0]
                        flag_lib[flag][multiple] += 1
            elif "not_labeled" in filters:
                skip = True
                try:
                    flag_lib["not_labeled, deleted"][0] += 1
                except KeyError:
                    flag_lib["not_labeled, deleted"] = [1,0]
            else :
                try:
                    flag_lib["not_labeled"][0] += 1
                except KeyError:
                    flag_lib["not_labeled"] = [1,0]
            if not skip:
                new_count += 1
        if not skip:
            #write line
            out.write(line+"\n")
    out.close()
    if progress:
        update_bash_line("writing new fastq file: "+str(i+1)+" / "+str(leng))
        print ("\tDone")
    return flag_lib,count,new_count

def flag_fastq(fastq,targets,progress,flags = []):#fastq file to be filtered, new file name, sequenses to be left out)
    leng = len(fastq)
    flag_lib = {}
    flag_lib["not_labeled"] = [0,0]
    start = True
    #if fastq.split(".")[-1] == "fastq":
    #   start = fastq[0].split(":")[0]
    #elif fastq.split(".")[-1] == "fasta":
    #    start = ">"
    if targets != []:
        for i in range(leng):
            line = fastq[i].strip()
            if start == True and line[0] == '@':
                start = start = fastq[0].split(":")[0]
            if i % 100 == 0 and progress:
                update_bash_line("labeling lines in fastq file: "+str(i)+" / "+str(leng))
            if line[:len(start)] == start:
                count = 0
                for x in range(len(targets)):
                    target = targets[x]
                    try:
                        addit = target[line[1:].split("|")[0].strip()]
                        if flags == []:
                            flag_star = addit[0]
                        else:
                            flag_star = flags[x]
                        if flag_star not in line.split():
                            line += " | " + flag_star +" - "+ addit[1]
                            count += 1
                        else:
                            None #when the flag is alreaddy present it is skipped
                    except KeyError:
                        None #this happens if the sequence was not alligned (so allot)
                if fastq[i] != line:
                    fastq[i] = line+'\n'
    if progress:
        update_bash_line("labeling lines in fastq file: "+str(i+1)+" / "+str(leng))
        print ("\tDone")
    return fastq

def write_report(report_file,report,file_number):
    file_log = {}
    with open(report_file,"r") as infile:
        for line in infile:
            title = line.strip().split('\t')[0]
            file_log[title] = line
    #create colmn names
    columns = file_log["file_name"].strip().split('\t')
    for key in report.keys():
        if key not in columns:
            columns.append(key)
    lines = list(file_log.keys())
    sorted(lines, reverse = True)
    #write file
    with open(report_file,"w") as infile:
        infile.write("\t".join(columns)+"\n")
        for line_key in lines:
            if line_key == "file_name":
                None
            elif line_key == file_number:
                line = file_log[line_key].strip()
                new_values = [str(report[x][0]+report[x][1]) for x in columns if x in report]
                infile.write(line+'\t'+'\t'.join(new_values)+'\n')
            else:
                infile.write(file_log[line_key])
    
def print_report(report,total,new_count,console_report):
    if console_report != "":
        with open(console_report,'w') as report_doc:
            report_doc.write("")
            report_doc.write('{0:^19}|{1:^19}|{2:^19}|{3:^19}'.format('flag','single','best of multi',"of total"))
            report_doc.write('-'*79)
            for key in report.keys():
                s = report[key][0]
                m = report[key][1]
                report_doc.write('{0:^19}|{1:^19}|{2:^19}|{3:^19}'.format(key,str(s),str(m),round(((s+m)/float(total)),2)))
            report_doc.write("")
            report_doc.write("number of seq: previous "+str(total)+", now "+str(new_count))
            report_doc.write("")
    print ("")
    print ('{0:^19}|{1:^19}|{2:^19}|{3:^19}'.format('flag','single','best of multi',"of total"))
    print ('-'*79)
    for key in report.keys():
        s = report[key][0]
        m = report[key][1]
        print ('{0:^19}|{1:^19}|{2:^19}|{3:^19}'.format(key,str(s),str(m),round(((s+m)/float(total)),2)))
    print ("")
    print("number of seq: previous "+str(total)+", now "+str(new_count))
    print ("")

def update_bash_line(string):
    stdout.write("\r"+string)
    stdout.flush()

def help():
    print ("script to flag sequenses of a fastq with sam file alignments" )
    print ("")
    print (" run as: python reads_filter.py -s [sam_file1,sam_file2] -i [fastq/fasta_file] -n [new_fasta/fastq]\n")
    print (" optional:\n-l [lable1,label2] :\tlabels to give to sam files\n-f [label1,label2] :\tlabels that will be removed\n-r [file_name] :\tmakes an report of removed reads\n-help or -h")
    print ("")

def main(arguments):
    print ("")
    if len(arguments) == 0:
        help()
    else:
        #standard options:
        sam_files = "in.sam"
        fastq = "in.fastq"
        output = "out.fastq"
        console_report = ""
        progress = True
        
        report_doc = ""
        flags = []
        filters = []
        #checking arguments for settings:
        for i in range(len(arguments)):
            arg = arguments[i]
            if arg == "-s": #sam file
                sam_files = arguments[i+1].strip().split(',')
            if arg == "-n": #new file
                output = arguments[i+1].strip()
            if arg == "-i": #in file
                fastq = arguments[i+1].strip()
            if arg == "-l": #label
                flags = arguments[i+1].strip().split(',')
            if arg == "-h" or arg == "-help":
                help()
            if arg == "-o": #console output writen in file 
                console_report = arguments[i+1].strip(',')
            if arg == "-p": #progress updates off
                progress = False
            if arg == "-f": #filter
                filters = [x.lower() for x in arguments[i+1].strip().split(",")]
            if arg == "-r": #write out removed reads (pipeline)
                report_doc = arguments[i+1].strip()
        #open fastq
        sample_number = fastq.split('/')[-1].split('_')[1]
        print ("reading "+fastq.split('/')[-1])
        fastq_open = open(fastq,"r").readlines()
        #label / filter
        if sam_files != "in.sam":
            targets = []
            for i in range(len(sam_files)):
                mapped_seq = read_sam(sam_files[i], progress)
                targets.append(mapped_seq)
            fastq = flag_fastq(fastq_open,targets,progress,flags)
        else:
            fastq = open(fastq,'r').readlines()
        #report / write
        (report,total,new_count) = stats_and_write(fastq,output,progress,filters)
        if report_doc != "":
            write_report(report_doc,report,sample_number)
        print_report(report,total,new_count,console_report)
        

if __name__ == '__main__':
    main(argv[1:])
