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

if not os.path.isfile('dgswe'):
  os.chdir(dgswe_work_dir)
  print "Compling dgswe"
  show_output('make clean')
  show_output('make dgswe')
  os.rename(dgswe_work_dir+'/'+'dgswe',example_dir+'/dgswe')

if not os.path.isfile('plot'):
  os.chdir(plot_work_dir)
  print "Compling plotting code"
  show_output('make clean')
  show_output('make plot')
  os.rename(plot_work_dir+'/'+'plot',example_dir+'/plot')

os.chdir(example_dir)

print "Running dgswe"
show_output('./dgswe')

print "Plotting results"
show_output('./plot')
