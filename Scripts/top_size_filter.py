"""
size filter
"""

from sys import (argv,stdout)


scrypt,infile,outfile,cutoff,report_file = argv

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
    sorted(lines)
    #write file
    with open(report_file,"w") as infile:
        infile.write("\t".join(columns)+"\n")
        for line_key in lines:
            if line_key == "file_name":
                None
            elif line_key == file_number:
                line = file_log[line_key].strip()
                new_values = [str(report[x]) for x in columns if x in report]
                infile.write(line+'\t'+'\t'.join(new_values)+'\n')
            else:
                infile.write(file_log[line_key])

count_filterd = 0
count = 0
print ("filtering "+infile.split('/')[-1])
with open(infile,'r') as open_infile:
    with open(outfile,'w') as open_outfile:
        write = True
        is_seq = True
        seq_info = ["",[],[]]
        start = ""
        for line in open_infile:
            if count % 100 == 0:
                stdout.write("\r"+"sequences checked: "+str(count))
                stdout.flush()
            if start == "" and line[0] == "@":
                start = line.split(':')[0]
            if line[:len(start)] == start:
                count += 1
                if seq_info[0] != "":
                    if len(''.join(seq_info[1])) <= int(cutoff):
                        open_outfile.write('@'+seq_info[0])
                        open_outfile.write('\n'.join(seq_info[1])+'\n')
                        open_outfile.write('+'+seq_info[0])
                        open_outfile.write('\n'.join(seq_info[2])+'\n')
                    else:
                        count_filterd += 1
                seq_info = [line[1:],[],[]]
                is_seq = True
            elif line[:len(start)] == '+'+start[1:]:
                is_seq = False
            else:
                if is_seq:
                    seq_info[1] += [line.strip()]
                else:
                    seq_info[2] += [line.strip()]
        #printing the last one
        if len(''.join(seq_info[1])) <= int(cutoff):
            open_outfile.write('@'+seq_info[0])
            open_outfile.write('\n'.join(seq_info[1])+'\n')
            open_outfile.write('+'+seq_info[0])
            open_outfile.write('\n'.join(seq_info[2])+'\n')
    #and done
    stdout.write("\r"+"sequences checked: "+str(count))
    stdout.flush()
    print ("\tdone")
    write_report(report_file,{"second_top_filter":count_filterd},infile.split('/')[-1].split('_')[1])
    print ("number filtered: "+str(count_filterd))

