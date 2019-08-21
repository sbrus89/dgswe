import subprocess
import os

def show_output(command):
  output = subprocess.Popen(command,stdout=subprocess.PIPE,shell=True)
  while output.poll() is None:
    l = output.stdout.readline()
    print l.rstrip('\n')

plot_work_dir = os.path.abspath('../../plot/work')
dgswe_work_dir = os.path.abspath('../../work')
example_dir = os.path.abspath('.')
nprocessors = '192'

if not os.path.isfile('dgswe_mpi'):
  os.chdir(dgswe_work_dir)
  print "Compling dgswe"
  show_output('make clean')
  show_output('make metis')
  show_output('make dgswe OPTS=mpi')
  show_output('make dgprep')
  show_output('make dgpost')
  os.rename(dgswe_work_dir+'/'+'dgswe_mpi',example_dir+'/dgswe_mpi')
  os.rename(dgswe_work_dir+'/'+'dgprep',example_dir+'/dgprep')
  os.rename(dgswe_work_dir+'/'+'dgpost',example_dir+'/dgpost')

if not os.path.isfile('plot'):
  os.chdir(plot_work_dir)
  print "Compling plotting code"
  show_output('make clean')
  show_output('make plot')
  os.rename(plot_work_dir+'/'+'plot',example_dir+'/plot')

os.chdir(example_dir)

print "Running dgswe"
f = open('np.in','w')
f.write(nprocessors)
f.close()
show_output('./dgprep < np.in')
show_output('mpirun -n '+nprocessors+' ./dgswe_mpi')
show_output('./dgpost')


print "Plotting results"
show_output('./plot')
