import subprocess

subprocess.call(['rm *.sol'],shell=True)
subprocess.call(['rm *.jpg'],shell=True)
subprocess.call(['rm *.out'],shell=True)
subprocess.call(['rm *.d'],shell=True)
subprocess.call(['rm *.ps'],shell=True)
subprocess.call(['rm fort.*'],shell=True)
subprocess.call(['rm -r PE*'],shell=True)
subprocess.call(['rm np.in'],shell=True)

