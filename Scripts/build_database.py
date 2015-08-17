
from sys import argv
import os

script,infile,outfile_dummy = argv
outfile = '.'.join(outfile_dummy.split('.')[:-1])
bash_command = "bowtie-build -f "+infile+" "+outfile
os.system(bash_command)
with open(outfile_dummy,'w') as write_file:
    write_file.write("dummy file")
