import subprocess
import sys
import os

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

f = open("../src/version.F90","r+")
lines = f.read().splitlines()

#print len(sys.argv)
#print sys.argv

git_branch = '"' + sys.argv[1] + '"' 
gitSHA = '"' + sys.argv[2] + '"'
compiler_version = '"' + sys.argv[3] + '"'
compiler_flags = '"' + sys.argv[4] + '"'
modified_files = '"' + sys.argv[5] + '"'
compile_date = '"' + sys.argv[6] + '"'
host = '"' + sys.argv[7] + '"'

lines[24] = "      gitBranch = " + git_branch 
lines[25] = "      gitSHA = " + gitSHA 
lines[26] = "      compiler_version = " + compiler_version
lines[27] = "      compiler_flags =  " + compiler_flags
#lines[28] = "     modified_files = " + modified_date
lines[29] = "      compile_date =  " + compile_date
lines[30] = "      host =  " + host

f.seek(0)
f.truncate()
for line in lines:
  f.write(line+"\n")
#  print line
f.close()
