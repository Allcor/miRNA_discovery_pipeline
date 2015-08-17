#! /usr/bin/python

"""
A python script that reads filtered fastq to take out the
sequences labeled as miRNA. It also counts them and assigns
names to them according to the alignment to miRBase by
reading the respective sorted .bam file (as .txt).
To do this, the script takes the sequences from the fastq
to which the "miRNA" label has been added. It puts these in
dictionary that is like: {Sequence: # times it occurs}.
These sequences are then looked up in the .txt file to add
the name of the reference sequence to the value of the dict-
ionary. In the end all information is sorted accoring to the
number reads and written to a fasta file.
"""

from sys import argv

def read_fastq(fastq):
    """
    input = fastq with labels
    output = dictionary with miRNA sequences, number, name,
    and alignment score
    
    A function that reads the fastq file, looks at 'miRNA'
    labels. If there are other labels as well, it takes the
    best possible score among them. Put the sequences that
    belong to miRNAs in a dictionary with their number (of
    occurence), name and the miRNA alignment score.    
    """
    sequence_dict = {}
    with open(fastq, 'r') as fastqfile:
        for line in fastqfile:
            line = line.strip('\n')
            ## read the file and strip newlines from it
            
            sequence = ''
            name = ''
            ## both set initially to empty string
            miRNAscore = 50
            rRNAscore = 50
            tRNAscore = 50
            
            ## Set scores for RNA species so that it is
            ## initially too high to be considered 'good'.
            ## This way, it won't interfere with the comparison
            
            if line[0] == '@':
                name = line
                ## if the line starts with @ it is the name
                
            if 'miRNA' in name:
            ## if 'miRNA' is in the name, it has been labeled
            ## as miRNA
                name = line.split('|')[1].strip()
                if 'rRNA' in name or 'tRNA' in name:
                    nameparts = name.split('|')
                    ## but it may have also been labeled as
                    ## other RNA species...
                    
                    for part in nameparts:
                        if 'miRNA' in part:
                            miRNAscore = int(part[8:])
                        elif 'rRNA' in part:
                            rRNAscore = int(part[7:])
                        elif 'tRNA' in part:
                            tRNAscore = int(part[7:])
                    ## check each of the alignment scores
                    if min(miRNAscore, rRNAscore, tRNAscore) == miRNAscore:
                        name = "miRNA - " + str(miRNAscore)
                        
                ## if the best score is miRNA, then the
                ## name is altered te reflect that
                
                    else:
                        name = ''
                        ## otherwise, empty the name (it is
                        ## then rRNA or tRNA)
                        
                nextline = fastqfile.readline()
                sequence = nextline.strip('\n')
                
                if len(sequence) > 0 and len(name) > 0:
                        
                    if sequence in sequence_dict:
                        sequence_dict[sequence][0] += 1
                    else:
                        sequence_dict[sequence] = [1, name]
                    
    print("Done with reading fastq!")
    return sequence_dict

def check_sam(SAM, sequence_dict):
    """
    input = sorted bamfile (as .txt) with alignments 
    (seq reads -> miRBase),
    and sequence_dict from the function "read_fastq" above
    output = list of all information that we need as:
    [number of reads, 'reference name', 'sequence']
    (list of lists)
    
    This function takes the sorted bam file (.txt) and 
    the dictionary of sequences and alignment scores to
    append the right name to each element in the dictionary,
    pop the dictionary values, discard sequences that occur
    only once and put the information in a list.
    """
    writelist = []
    with open(SAM, 'r') as samfile:
        for line in samfile:
        # read the file
            name = line.split('\t')[2]
            # write down the name
            sequence = line.split('\t')[9]
            # and the sequence
            for key in sequence_dict:
                if key == sequence:
                # if the key is the sequence you found in the sam
                    sequence_dict[key].append(name)
                    # add the name to it
            if sequence in sequence_dict:
                if len(sequence_dict[sequence]) > 2:
                    write = sequence_dict.pop(sequence)
                    if int(write[0]) > 1:
                    ## only take sequences that occur more than once
                    ## (you get loads of crap if you keep 'em all)
                        writelist.append([write[0], write[2], sequence])
                    
    print("The list is now ready to be sorted and written!")
    return writelist

def collapse(writelist):
    writelist.sort(key=lambda x: x[1])
    writelist_collapsed = []
    same_ref = []
    same_name = ""
    for seq in writelist:
        if same_name != seq[1]:
            if same_name != "":
                same_ref.sort(key=lambda x: x[0], reverse = True)
                sequence = same_ref[0][2]
                count = 0
                for i in same_ref:
                    count += i[0]
                writelist_collapsed.append([count, same_name, sequence])
            same_ref = [seq]
            same_name = seq[1]
        else:
            same_ref.append(seq)
    #writing last line
    same_ref.sort(key=lambda x: x[0], reverse = True)
    sequence = same_ref[0][2]
    count = 0
    for i in same_ref:
        count += i[0]
    writelist_collapsed.append([count, same_name, sequence])
    #end
    return writelist_collapsed
           
def write_to_file(writelist, OUTPUT):
    """
    input = list of all the information to be written as:
    [nr of reads, 'reference sequence name', 'sequence']
    and the file to which this information is to be written
    output = fasta file with sequences sorted by occurence
    
    This function uses the writelist made by check_sam (above)
    , sorts it and writes it to the fasta file specified at
    the top.
    """
    writelist.sort(key=lambda x: x[0], reverse = True)
    ## sort the list by number of reads, descending
    ## (highest first)
    number_of_miRNAs = 0
    with open(OUTPUT, 'w') as writefile:
        for lists in writelist:
            number_of_miRNAs += lists[0]
            writefile.write('>' + lists[1] + "_x" + str(lists[0]) + '\n' + lists[2] + '\n')
    print(number_of_miRNAs, "miRNAs have been identified and written to .fasta.")
        

def main(arguments):
    fastq = arguments[0]
    sam = arguments[1]
    output = arguments[2]
    write_to_file(collapse(check_sam(sam, read_fastq(fastq))), output)
    
    
if __name__ == '__main__':
    main(argv[1:])
