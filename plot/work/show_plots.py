import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib import interactive
import numpy as np
import os
import curses

interactive(True)

frmt = 'jpg'
sol = 'zeta'

direc_ls = [
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_tri_x64/p3/ctp3/hbp3/plots/upper/',
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_tri/p3/ctp3/hbp3/plots/upper/',
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_tri/p3/ctp1/hbp1/plots/upper/',
            
            '/home/sbrus/data-drive/galveston_spline_flux/galveston_tri_x64/p3/ctp3/hbp3/plots/inlet/',
            '/home/sbrus/data-drive/galveston_spline_flux/galveston_tri/p3/ctp3/hbp3/plots/inlet/',
            '/home/sbrus/data-drive/galveston_spline_flux/galveston_tri/p3/ctp1/hbp1/plots/inlet/',
            
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_tri_x64/p3/ctp3/hbp3/plots/left_inlet/',
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_tri/p3/ctp3/hbp3/plots/left_inlet/',
            #'/home/sbrus/data-drive/galveston_spline_flux/galveston_tri/p3/ctp1/hbp1/plots/left_inlet/'            
            
            ]         

         
keep_going = True  
snap = 1
while keep_going:          

  snap_str = "%04d" % snap
  
  for i,direc in enumerate(direc_ls):         

    im_file = direc + sol + '_' + snap_str + '.' + frmt 
    if not os.path.isfile(im_file):
      keep_going = False
      break
  
    fig = plt.figure(i)
    fig.canvas.set_window_title(direc)
    img=mpimg.imread(im_file)      
    plt.imshow(img)
    plt.tight_layout()
    plt.axis('off')


  inp = raw_input("Next image: >, Previous image: <  ")   
  while inp != ',' and inp != '.':
    inp = raw_input("Next image: >, Previous image: <  ") 
    
  if inp == '.':
    snap = snap + 1
  elif inp == ',' and snap > 1:
    snap = snap - 1

