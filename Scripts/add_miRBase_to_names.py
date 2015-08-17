#! /usr/bin/python

'''
Make a list of miRNA names and sequences, to be able to
allocate the right names to the right miRNAs. Takes the names
by Barozai et al. (2013) as starting point and fills the list
with miRNAs found with miRBase.
'''

NAMES_LIST = "/zfs/datastore0/group_root/GreenStudentLab/001-miRNA_discovery/002-Carrot/pipeline/Barozai_miRNAs.txt"
MIRBASE_MIRNAS = "/zfs/datastore0/group_root/GreenStudentLab/001-miRNA_discovery/002-Carrot/pipeline/Results/miRBase_found-all.fasta"

NAMES_DICT = {}
with open(NAMES_LIST, 'r') as names:
# initialise dictionary with names from Barozai
    for line in names:
        line = line.strip('\n')
        name = line.split()[0]
        sequence = line.split()[1]
        NAMES_DICT[sequence] = name
        
with open(MIRBASE_MIRNAS, 'r') as new_names:
# open the miRBase_found-all fasta file
    for line in new_names:
        if line[0] == '>':
            line = line.strip('\n')[8:].split('_')[0]
            # cut out the name part if the line starts with
            # '>'
            name = ''
            for character in line:
                if not character.isalpha():
                    name += character
                else:
                    break
            # use only the numbers of the name
            name = "dca-MIR" + name
            # and append those numbers to "dca-MIR"
        else:
            sequence = line.strip('\n').replace('T', 'U')
        # if the line does not start with a '>', it is the
        # sequence
        if sequence in NAMES_DICT:
            pass
        # if the sequence already exists in the dictionary
        # it is already there and does not need to be added
        # again
        else:
        # if the sequence is NOT in the dictionary, it needs
        # to be added and the names may have to be adjusted
            NAMES_DICT[sequence] = name


sort_list = []
for value in NAMES_DICT.values():
    sort_list.append(value)
    sort_list.sort()
    
print(sort_list)
