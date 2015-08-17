#! /usr/bin/python

from sys import argv
from tex_table import *

script, ath, nr = argv

if ath == 'y':
    inputfile = "/zfs/datastore0/group_root/GreenStudentLab/001-miRNA_discovery/002-Carrot/Results/Drylab-data/Bowtiemapping/Sample" + nr + "_miRDP_mature-miRBase_ath.sam"
    outputfile = "/zfs/datastore0/group_root/GreenStudentLab/001-miRNA_discovery/002-Carrot/Results/Drylab-data/Tables/Sample" + nr + "_miRNA-table_mature_ath.csv"
elif ath == 'n':
    inputfile = "/zfs/datastore0/group_root/GreenStudentLab/001-miRNA_discovery/002-Carrot/Results/Drylab-data/Bowtiemapping/Sample" + nr + "_miRDP_mature-miRBase.sam"
    outputfile = "/zfs/datastore0/group_root/GreenStudentLab/001-miRNA_discovery/002-Carrot/Results/Drylab-data/Tables/Sample" + nr + "_miRNA-table_mature.csv"
elif ath == 'nim':
    inputfile = "/zfs/datastore0/group_root/GreenStudentLab/001-miRNA_discovery/002-Carrot/Results/Drylab-data/Bowtiemapping/Sample" + nr + "_miRDP_mature-miRBase+unal.sam"
    outputfile = "/zfs/datastore0/group_root/GreenStudentLab/001-miRNA_discovery/002-Carrot/Results/Drylab-data/Tables/Sample" + nr + "_mature_notinmiRBase.csv"


hairpins = "/zfs/datastore0/group_root/GreenStudentLab/001-miRNA_discovery/002-Carrot/Results/Drylab-data/miRDeep/Sample"+ nr +"_filter_P_prediction"
predictions = "/zfs/datastore0/group_root/GreenStudentLab/001-miRNA_discovery/002-Carrot/Results/Drylab-data/miRDeep/Sample" + nr + "_predictions"

# inputfile = bowtie2 output: Sample#_miRDP_mature-miRBase*.sam
# outfile = table to be made (.csv?)
# hairpins = Sample#_filter_P_predictions
# and predictions = Sample#_predictions

def make_table(inputfile, outputfile):
    miRNA_family_members = {}
    hairpin_dict = {}
    with open(hairpins, 'r') as hairpinfile:
        for line in hairpinfile:
            line = line.split()
            name = line[2]
            PL = line[7]
            hairpin_dict[name] = PL
            
    with open(inputfile, 'r') as readfile:
        with open(outputfile, 'w') as writefile:
        #write a header to the outfile
            writefile.write("Source miRNA\tmiRNA\tPrec Len\tMin Free Energy\tSequence\tMism\tLength\tMat Seq Arm\tNo. Reads\n")
            for line in readfile:
            # read the whole file, line by line
                line = line.strip('\n')
                # and take out "extra lines"
                if line[0] != '@':
                # ignore headers
                    line = line.split('\t')
                    # split each line by tabs
                    name,source = columns1_2(line)
                    # name = new miRNA name
                    #source = 'source miRNA name'
                    
                    PL = len(hairpin_dict[line[0]])
                    # PL = precursor length (= column 3)
                    
                    mfe = column4(line, predictions)[0]
                    # mfe = minimum free energy of pre-miRNA molecule,
                    # as calculated by RNAfold
                    
                    sequence, mismatches, mat_length = columns5_6_7(line) #line[9]
                    # sequence = mature sequence
                    # mismatches = mismatches in mature sequence, compared to source
                    # mat_length = length of mature miRNA
                                                    
                    MSA = column9(line, predictions)
                    # mature sequence arm (which arm of the hairpin the miRNA is on
                    
                    nr_of_reads = column4(line, predictions)[1]
                    nr_of_reads = int(nr_of_reads.split('x')[1])
                    
                    if source not in miRNA_family_members.keys():
                        miRNA_family_members[source] = [name, str(PL), mfe, sequence, str(mismatches), str(mat_length), MSA, str(nr_of_reads),hairpin_dict[line[0]]]
            ## This creates the dict that contains all info as:
            ## {'source' : ['new_name', 'PL', 'mfe', 'seq', 'mm', 'mat_len', 'msa', 'nr_reads', 'hairpin_seq']}
                    
            #### WRITE WITHOUT DUPLICATES, SORTED
            sortlist = []
            for key in miRNA_family_members.keys():
                sortlist.append(key)
            
            sortlist.sort()
            
            for keys in sortlist:                
                writefile.write(keys +'\t' + '\t'.join(miRNA_family_members[keys]) + '\n')
                Tex_graph.row([keys]+miRNA_family_members[keys])
            Tex_graph.unique_names()
            Tex_graph.write("test.tex")

def columns1_2(line):
    source = line[2]
    memberID = source.split('-')[1]
    memberID = ''.join(i for i in memberID if i.isdigit())
    name = "dca-miR" + memberID

    return name, source
    
def column4(line, predictions):
    name = line[0]
    write = False
    with open(predictions, 'r') as pred_file:
        for pred_line in pred_file:
            pred_line = pred_line.strip('\n')
            if name in pred_line:
                write = True
            if write:
                if 'mfe' in pred_line:
                    return pred_line.split()[1], name
                
def columns5_6_7(line):
    sequence = line[9]
    mat_length = len(sequence)
    new_sequence = ''
    mismatches = 0
    
    for element in line:
        if element[:4] == 'MD:Z':
            MDZ = element[5:]
            # find the MD:Z value and make the variable
            
    MDZc = ''
    # take a copy of MD:Z, to easily fetch remaining numbers
    if len(MDZ) > 2:
        for character in MDZ:
            if character not in 'ACGT^':
                MDZc += character
            if character in 'ACGT^':
                if MDZc == '':
                    MDZc = 0
                mismatches += 1
                lenseq = len(new_sequence)
#                if character == '^':
#                    new_sequence += sequence[lenseq:lenseq+int(MDZc)] + '-'
#                else:
#                    new_sequence += sequence[lenseq:lenseq+int(MDZc)] + character.lower()
#                MDZc = ''

### NOTE: this way you are recreating the reference sequence!!
### If you want to keep the carrot sequence, you should not replace
### the nucleotides in the MDZ.

                if character == '^':
                ### this is only 1 sequence... (bna-miR169l)
                    new_sequence += sequence[lenseq:lenseq+int(MDZc)+1]
                else:
                    new_sequence += sequence[lenseq:lenseq+int(MDZc)] + sequence[lenseq+int(MDZc)].lower()
                MDZc = ''
                
### This should keep the carrot sequence, but change the mismatches
### to lowercase letters...

        lenseq = len(new_sequence)
        new_sequence = new_sequence + sequence[lenseq:]
        new_sequence = new_sequence.replace('T', 'U').replace('t', 'u')
        
        ## also change the DNA sequence back into RNA
        ## (makes sence for mi*RNA*s, right?)
        
    else:
        new_sequence = sequence.replace('T', 'U')
        
    return new_sequence, mismatches, mat_length
    
def column9(line, predections):
    name = line[0]
    seq_arm = ''
    
    with open(predictions, 'r') as pred_file:
        for pred_line in pred_file:
            pred_line = pred_line.strip('\n')
            if pred_line != '':
                if name in pred_line:
                    if seq_arm == "first":
                        return "5'"
                    elif seq_arm == "second":
                        return "3'"
                    else:
                        return seq_arm
                if pred_line.split()[0] == 'mature_arm':
                    seq_arm = pred_line.split()[1]
                    
#################################
## Write similar functions for miRNAs that are not in miRBase.
## Notice that in that case, you lack:
## 1. Source miRNA
## 2. Mismatches
##
## Therefore, you can also not make the names in the same way.
## So, use the name "dca-MIR???" for now.
## The sequence does not have to be edited.
## Also, sequences now have to be filtered by the "4" in the
## .sam file.
##
## The columns are then:
## 1 - name = "dca-MIR???"
## 2 - precursor_length
## 3 - minimum free energy
## 4 - sequence
## 5 - mature sequence length
## 6 - mature sequence arm
##
#################################

def col1(line):
    return "dca-miR???"

def col3(line, predictions):
    name = line[0]
    write = False
    with open(predictions, 'r') as pred_file:
        for pred_line in pred_file:
            pred_line = pred_line.strip('\n')
            if name in pred_line:
                write = True
            if write:
                if 'mfe' in pred_line:
                    return pred_line.split()[1]

def col4_5(line):
    sequence = line[9]
    matseqlen = len(sequence)
    sequence = sequence.replace('T', 'U')
    return sequence, matseqlen

def col6(line, predictions):
    name = line[0]
    seq_arm = ''
    
    with open(predictions, 'r') as pred_file:
        for pred_line in pred_file:
            pred_line = pred_line.strip('\n')
            if pred_line != '':
                if name in pred_line:
                    if seq_arm == "first":
                        return "5'"
                    elif seq_arm == "second":
                        return "3'"
                    else:
                        return seq_arm
                if pred_line.split()[0] == 'mature_arm':
                    seq_arm = pred_line.split()[1]
                    
def make_table_not_in_miRBase(inputfile, outputfile):
    sortdict = {}
    hairpin_dict = {}
    with open(hairpins, 'r') as hairpinfile:
        for line in hairpinfile:
            line = line.split()
            name = line[2]
            PL = len(line[7])
            hairpin_dict[name] = PL
            
    with open(inputfile, 'r') as readfile:
        with open(outputfile, 'w') as writefile:
        #write a header to the outfile
            writefile.write("miRNA\tPrec Len\tMin Free Energy\tSequence\tLength\tMat Seq Arm\tNo. Reads\n")
            refID = 0
            for line in readfile:
                line = line.strip('\n').split('\t')
                if line[1] == '4':
                    refID += 1
                    name = col1(line)
                    prec_len = hairpin_dict[line[0]]
                    mfe = col3(line, predictions)
                    sequence, length = col4_5(line)
                    matseqarm = col6(line, predictions)
                    nr_of_reads = line[0].split('x')[1]
                    
                    ## write all to a dictionary, to sort
                    
                    sortdict[refID] = [name, prec_len, mfe, sequence, matseqarm, nr_of_reads]
                    
                    ## print for testing purposes
#                    print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (name, prec_len, mfe, sequence, length, matseqarm, nr_of_reads)

                    ## write for actual use - old style
                    #writefile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (name, prec_len, mfe, sequence, length, matseqarm, nr_of_reads))
## now write the sorted list!
            sortlist = []
            for key in sortdict.keys():
                sortlist.append((key, int(sortdict[key][5])))
            sortlist.sort(key=lambda tup: tup[1])
            sortlist = sortlist[::-1]
            for keys in sortlist:
                writefile.write('\t'.join(str(x) for x in sortdict[keys[0]]) + '\n')

### Options to make different tables, depending on if the reads
### have been aligned to miRBase or not.
### if ath = 'y' or 'n', run the table for the miRNAs in miRBase.
### if it is 'nim' (for Not In MiRBase), make the other table.

Tex_graph = miR_table()

if ath == 'y' or ath == 'n':
    make_table(inputfile, outputfile)
elif ath == 'nim':
    make_table_not_in_miRBase(inputfile, outputfile)
