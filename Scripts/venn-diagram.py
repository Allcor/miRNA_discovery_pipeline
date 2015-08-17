#!/usr/bin/python

from sys import argv

script,barozai_file,mirbase_file,miRdeep_file,outfile = argv
#barozai_file = "/zfs/datastore0/group_root/GreenStudentLab/001-miRNA_discovery/002-Carrot/Results/Drylab-data/PakiBase(copies)/pakibase-venn.fasta"
#mirbase_file = "./Results/miRBase_found-all.fasta"
#miRdeep_file = "./Report/Tables/miRDP_found-all.fasta"


dictio = {}
    
def read_fasta(filename,tag):
    with open(filename,'r') as infile:
        for line in infile:
            if line[0] != '>':
                if line.strip() in dictio:
                    dictio[line.strip()] += [tag]
                else:
                    dictio[line.strip()] = [tag]
                    
read_fasta(barozai_file,"baro")
read_fasta(mirbase_file,"base")
read_fasta(miRdeep_file,"deep")

in_all = 0
baro_base = 0
deep_base = 0
baro_deep = 0
baro = 0
deep = 0
base = 0

for seq in dictio:
    venns = dictio[seq]
    if "baro" in venns:
        if "base" in venns:
            if "deep" in venns:
                in_all += 1
            else:
                baro_base += 1
        elif "deep" in venns:
            baro_deep += 1
        else:
            baro += 1
    elif "base" in venns:
        if "deep" in venns:
            deep_base += 1
        else:
            base += 1
    elif "deep" in venns:
        deep += 1
        
with open(outfile,'w') as out_file:
    outfile.write("barozai;mirdeep;mirbase;Counts")
    outfile.write("0;0;0;"+str(baro+baro_deep+baro_base+deep+deep_base+base+in_all))
    outfile.write("0;0;1;"+str(base))
    outfile.write("0;1;0;"+str(deep))
    outfile.write("0;1;1;"+str(deep_base))
    outfile.write("1;0;0;"+str(baro))
    outfile.write("1;0;1;"+str(baro_base))
    outfile.write("1;1;0;"+str(baro_deep))
    outfile.write("1;1;1;"+str(in_all))
