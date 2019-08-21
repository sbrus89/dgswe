import subprocess
import os

plot_work_dir = os.path.abspath('../../plot/work')
dgswe_work_dir = os.path.abspath('../../work')
example_dir = os.path.abspath('.')

if not os.path.isfile('dgswe'):
  os.chdir(dgswe_work_dir)
  print "Compling dgswe"
  subprocess.call(['make','clean'])
  output = subprocess.Popen(['make','dgswe'],stdout=subprocess.PIPE)
  while output.poll() is None:
    l = output.stdout.readline()
    print l.rstrip('\n')
  os.rename(dgswe_work_dir+'/'+'dgswe',example_dir+'/dgswe')

if not os.path.isfile('plot'):
  os.chdir(plot_work_dir)
  print "Compling plotting code"
  subprocess.call(['make','clean'])
  output = subprocess.Popen(['make','plot'],stdout=subprocess.PIPE)
  while output.poll() is None:
    l = output.stdout.readline()
    print l.rstrip('\n')
  os.rename(plot_work_dir+'/'+'plot',example_dir+'/plot')

os.chdir(example_dir)

print "Running dgswe"
output = subprocess.Popen(['./dgswe'],stdout=subprocess.PIPE)
while output.poll() is None:
  l = output.stdout.readline()
  print l.rstrip('\n')


print "Plotting results"
output = subprocess.Popen(['./plot'],stdout=subprocess.PIPE)
while output.poll() is None:
  l = output.stdout.readline()
  print l.rstrip('\n')
