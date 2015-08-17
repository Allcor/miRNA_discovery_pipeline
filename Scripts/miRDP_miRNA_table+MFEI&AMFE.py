#! /usr/bin/python

from sys import argv
from tex_tablePLUS import *

script, in_mirbase, inputfile, outputfile, hairpins, predictions = argv


# inputfile = bowtie2 output: Sample#_miRDP_mature-miRBase*.sam
# outfile = table to be made (.csv?)
# hairpins = Sample#_filter_P_prediction
# and predictions = Sample#_predictions

def make_table(inputfile, outputfile):
    miRNA_family_members = {}
    hairpin_dict = {}
    with open(hairpins, 'r') as hairpinfile:
        for line in hairpinfile:
            line = line.split()
            name = line[2]
            prec_seq = line[7]
            hairpin_dict[name] = prec_seq
            
    with open(inputfile, 'r') as readfile:
        with open(outputfile, 'w') as writefile:
        #write a header to the outfile
            writefile.write("Source miRNA\tmiRNA\tPrec Len\tMin Free Energy\tAMFE\tMFEI\tSequence\tMism\tLength\tMat Seq Arm\tNo. Reads\n")
            for line in readfile:
            # read the whole file, line by line
                line = line.strip('\n')
                # and take out "extra lines"
                if line[0] != '@':
                # ignore headers
                    line = line.split('\t')
                    # split each line by tabs
                    if line[1] == '4':
                        pass
                    else:
                        
                        name,source = columns1_2(line)
                        # name = new miRNA name
                        #source = 'source miRNA name'
                        
                        PL = len(hairpin_dict[line[0]])
                        # PL = precursor length (= column 3)
                        
                        mfe = str(abs(float(column4(line, predictions)[0])))
                        # mfe = minimum free energy of pre-miRNA molecule,
                        # as calculated by RNAfold
                        
                        amfe = (float(mfe) / PL) * 100
                        # Adjusted MFE, see Zhang, et al. 2006 (Evidence that miRNAs are different from other RNAs)
                        GChairpin = ((hairpin_dict[line[0]].count('C') + hairpin_dict[line[0]].count('G')) / PL * 100)
                        
                        mfei = (amfe / (GChairpin))
                        # MFE Index, also see Zhang, et al. 2006
                        
                        sequence, mismatches, mat_length = columns5_6_7(line) #line[9]
                        # sequence = mature sequence
                        # mismatches = mismatches in mature sequence, compared to source
                        # mat_length = length of mature miRNA
                                                        
                        MSA = column9(line, predictions)
                        # mature sequence arm (which arm of the hairpin the miRNA is on
                        
                        nr_of_reads = column4(line, predictions)[1]
                        nr_of_reads = int(nr_of_reads.split('x')[1])
                        
                        if source not in miRNA_family_members.keys():
                            miRNA_family_members[source] = [name, str(PL), mfe, format(amfe, '.2f'), format(mfei, '.2f'), sequence, str(mismatches), str(mat_length), MSA, str(nr_of_reads),hairpin_dict[line[0]]]
            ## This creates the dict that contains all info as:
            ## {'source' : ['new_name', 'PL', 'mfe', 'amfe', 'mfei' 'seq', 'mm', 'mat_len', 'msa', 'nr_reads', 'hairpin_seq']}
                    
            #### WRITE WITHOUT DUPLICATES, SORTED
            sortlist = []
            for key in miRNA_family_members.keys():
                sortlist.append(key)
            
            sortlist.sort()
            
            for keys in sortlist:                
                writefile.write(keys +'\t' + '\t'.join(miRNA_family_members[keys]) + '\n')
                #Tex_graph.row([keys]+miRNA_family_members[keys])
            #Tex_graph.unique_names()
            #Tex_graph.write("test.tex")

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
                    return [pred_line.split()[1], name]
                
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
            prec_seq = line[7]
            hairpin_dict[name] = prec_seq
            
    with open(inputfile, 'r') as readfile:
        with open(outputfile, 'w') as writefile:
        #write a header to the outfile
            writefile.write("miRNA\tPrec Len\tMin Free Energy\tAMFE\tMFEI\tSequence\tLength\tMat Seq Arm\tNo. Reads\n")
            refID = 0
            for line in readfile:
                line = line.strip('\n').split('\t')
                if line[1] == '4':
                    refID += 1
                    name = col1(line)
                    prec_len = len(hairpin_dict[line[0]])
                    mfe = str(abs(float(col3(line, predictions))))
                    amfe = format((float(mfe) / prec_len) * 100, '.2f')
                    GChairpin = ((hairpin_dict[line[0]].upper().count('C') + hairpin_dict[line[0]].upper().count('G')) / float(prec_len) * 100)
                    mfei = format((float(amfe) / (GChairpin)), '.2f')
                    sequence, length = col4_5(line)
                    matseqarm = col6(line, predictions)
                    nr_of_reads = line[0].split('x')[1]
                    
                    ## write all to a dictionary, to sort
                    
                    sortdict[refID] = [name, prec_len, mfe, amfe, mfei, sequence, length, matseqarm, nr_of_reads]
                    
                    ## print for testing purposes
#                    print "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (name, prec_len, mfe, amfe, mfei, sequence, length, matseqarm, nr_of_reads)

                    ## write for actual use - old style
                    #writefile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (name, prec_len, mfe, amfe, mfei, sequence, length, matseqarm, nr_of_reads))
## now write the sorted list!
            sortlist = []
            for key in sortdict.keys():
                sortlist.append((key, int(sortdict[key][8])))
            sortlist.sort(key=lambda tup: tup[1])
            sortlist = sortlist[::-1]
            for keys in sortlist:
                writefile.write('\t'.join(str(x) for x in sortdict[keys[0]]) + '\n')

### Options to make different tables, depending on if the reads
### have been aligned to miRBase or not.
### if ath = 'y' or 'n', run the table for the miRNAs in miRBase.
### if it is 'nim' (for Not In MiRBase), make the other table.

Tex_graph = miR_table()

if in_mirbase == 'y':
    make_table(inputfile, outputfile)
elif in_mirbase == 'n':
    make_table_not_in_miRBase(inputfile, outputfile)
else:
    print("Please enter a valid option - 'y' or 'n'.")
