
DATABASE_FILE = "new_miRNA_names.txt"

class Mir_row():
    def __init__(self,data):
        self.name = data[1]
        self.source = data[0]
        self.PL = data[2]
        self.MFE = data[3]
        self.AMFE = data[4]
        self.MFEI = data[5]
        self.MS = data[6]
        self.NM = data[7]
        self.ML = data[8]
        self.MSA = data[9]
        self.reads = data[10]
        self.pre = data[11]
    
    def write(self):
        line = [self.name,self.source,self.PL,self.MFE,self.AMFE, self.MFEI,self.MS,self.NM,self.ML,self.MSA]
        line = ' & '.join([str(x) for x in line])+'\\\\\n'
        return line

class miR_table():
    def __init__(self):
        self.data = []
        self.table_name = "Discovered miRNA"
        self.caption = "Name of the newly identified miRNA, closest miRNA found in miRBase, PL = Precursor miRNA Length, MFE = Minimum Free Energy given by RNAfold, NM = Number of mismatches relative to sequence of reference of column 2 (lower case is base found in \emph{Daucus carota} and not in reference), ML = Mature sequence Length, MSA Mature Sequence Arm."

    def row(self,columns):
        self.data.append(Mir_row(columns))
    
    def unique_names(self):
        #logging previous miRNA
        mi_file = open(DATABASE_FILE,'a+')
        mi_dic = {}
        mi_list = []
        for line in mi_file:
            line = line.strip().split('\t')
            name = line[0]
            precursor = line[1]
            frame = line[2]
            mi_dic[(precursor,frame)] = name
        mi_file.close()
        for row in self.data:
            try:
                row.name = mi_dic[(row.pre,row.MSA)]
            except KeyError:
                mi_list.append(row)
        #creating names
        names = [((mi_dic[key].split('-')[0]+'-'+mi_dic[key].split('-')[1])[:-1],mi_dic[key].split('-')[1][-1],key) for key in mi_dic.keys()]
        names.sort(key=lambda x: x[0]+x[1])
        if names == None:
            names = []
        for row in mi_list:
            count = 0
            x = 0
            for i in range(len(names)):
                if names[i][0] == row.name:
                    x = i
                    while names[x][0] == row.name:
                        if row.pre == names[x][2][0]:
                            if names[x][2][1] != row.MSA: #just making sure
                                prime = filter(lambda x: x.isdigit(), row.MSA)
                                print('double trouble! '+row.name+names[x][1]+'-'+prime+'p')
                                mi_dic[(row.pre,row.MSA)] = row.name+names[x][1]+'-'+prime+'p'
                                prime = filter(lambda x: x.isdigit(), names[x][2][1])
                                mi_dic[(names[x][2][0],names[x][2][1])] = mi_dic[(names[x][2][0],names[x][2][1])]+'-'+prime+'p'
                                names.insert(x,(row.name,names[x][1],row.MSA))
                                break
                            else:
                                print("IMPOSEBRU!!!")
                        elif names[x][2][0] != names[x+1][2][0]:#next precursor not the same
                            count += 1
                        else:
                            print (names[i])#should not happen
                        x += 1
                    break
            new_name = row.name+"abcdefghijklmno"[count]#gives error if more, not sure if they are suposed to go past 15
            print (new_name+' added')
            mi_dic[(row.pre,row.MSA)] = new_name
            names.insert(x,(new_name[:-1],new_name[-1],(row.pre,row.MSA)))
        #writing file again
        line_list = []
        for key in mi_dic.keys():
            line_list.append(mi_dic[key]+'\t'+key[0]+'\t'+key[1])
        line_list.sort()
        mi_file = open(DATABASE_FILE,'w')
        for line in line_list:
            mi_file.write(line+'\n')
        mi_file.close()
        #changing names in row class
        for x in self.data:
            x.name = mi_dic[(x.pre,x.MSA)]
        
    def edit(self,y,x,value):
        if x == name:
            self.data[y].name = value
        elif x == source:
            self.data[y].source = value
        elif x == pl:
            self.data[y].PL = value
        elif x == mfe:
            self.data[y].MFE = value
        elif x == ms:
            self.data[y].MS = value
        elif x == nm:
            self.data[y].NM = value
        elif x == ml:
            self.data[y].ML = value
        elif x == msa:
            self.data[y].MSA = value

    def write(self,filename):
        tex_file = open("./Scripts/latex/"+filename,"w")
        tex_file.write("\documentclass[preview=true,varwidth=\maxdimen]{standalone}\n")
        tex_file.write("\\usepackage{booktabs}\n")
        tex_file.write("\\begin{document}\n\\begin{table}\n\centering\n")
        tex_file.write("\caption{"+self.table_name+"}")
        tex_file.write("\\begin{tabular}{c c c c c c c c}\n\\\\\\toprule\n\\textbf{miRNAs} & \\textbf{Source miRNAs} & \\textbf{PL} & \\textbf{MFE} & \\textbf{Mature Sequence} & \\textbf{NM} & \\textbf{ML} & \\textbf{MSA}\\\\\n\midrule\n")
        for line in self.data:
            tex_file.write(line.write())
        tex_file.write("\\bottomrule\n\end{tabular}\n")
        tex_file.write("\label{tab:template}\n\end{table}\n\\begin{minipage}{16cm}\n")
        tex_file.write(self.caption)
        tex_file.write("\end{minipage}\n\end{document}")
        tex_file.close()
