import os
import sys
import shutil
import glob
import subprocess

def write_file(name,content):    
      
      
  if name != '':
    f = open(name,'w')   
    for line in content:
      value = line['value']
      comment = line['comment']    
      spaces = 101 - len(value)    
      f.write(value + spaces*' ' + comment)
    f.close()   
    
    

plot_exe_direc = '/home/sbrus/Codes/dgswe/plot/work/'

#results_direc = '/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri/p3/ctp1/hbp1/'
results_direc = '/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri/p3/ctp3/hbp3/'
#results_direc = '/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri/adcirc/ESL0/'
#results_direc = '/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri_x64/adcirc/ESL0/'

plotting_direc = results_direc+'plots/'
#plotting_direc = results_direc+'plots/dissertation/'
#plotting_direc = results_direc+'plots/presentation/'

zoom_list = [
             ['2.95e5,3.12e5,3.283e6,3.298e6','upper'],  
             ['3.2e5,3.4e5,3.24e6,3.26e6','inlet'],
             ['2.84e5,3.00e5,3.21e6,3.233e6','left inlet']
             ]
             
plot_order = '20,20,1'             
#plot_order = '3,3,1'              
             
if not os.path.exists(plotting_direc):
  os.makedirs(plotting_direc)
 
adcirc = False
if 'adcirc' in results_direc:
  adcirc = True
             
shutil.copy(plot_exe_direc+'plot',plotting_direc+'plot')
shutil.copy(results_direc+'dgswe.inp',plotting_direc)

if adcirc:
  subprocess.Popen('ln -s '+results_direc+'fort.63 '+plotting_direc, shell=True, stdout=subprocess.PIPE).communicate()[0] 
  subprocess.Popen('ln -s '+results_direc+'fort.64 '+plotting_direc, shell=True, stdout=subprocess.PIPE).communicate()[0]  
else:  
  subprocess.Popen('ln -s '+results_direc+'Z.sol '+plotting_direc, shell=True, stdout=subprocess.PIPE).communicate()[0] 
  subprocess.Popen('ln -s '+results_direc+'Qx.sol '+plotting_direc, shell=True, stdout=subprocess.PIPE).communicate()[0] 
  subprocess.Popen('ln -s '+results_direc+'Qy.sol '+plotting_direc, shell=True, stdout=subprocess.PIPE).communicate()[0] 
  subprocess.Popen('ln -s '+results_direc+'hb.sol '+plotting_direc, shell=True, stdout=subprocess.PIPE).communicate()[0] 

 

for zoom in zoom_list :  
  
  zoom_box  = zoom[0]
  zoom_name = zoom[1]
  zoom_name_ = "_".join(zoom_name.split())
  
  save_direc = plotting_direc+zoom_name_+'/'
  if not os.path.exists(save_direc):
    os.makedirs(save_direc)
  
  
  base_cscale_direc = '/home/sbrus/data-drive/galveston_spline_flux_fix/galveston_tri_x64/p3/ctp3/hbp3/plots/'+zoom_name_+'/'
  shutil.copy(base_cscale_direc+'vel.cscale',plotting_direc)
  shutil.copy(base_cscale_direc+'zeta.cscale',plotting_direc)

  plot_input = [{'value':'./'                                             , 'comment':'!'+'dgswe.inp file path \n'},
                {'value':'/scratch365/sbrus/,/home/sbrus/data-drive/'     , 'comment':'!'+'dgswe.inp path substitution \n'},
                {'value':'1,2,1'                                          , 'comment':'!'+'zeta plot options \n'},
                {'value':'1,2,1'                                          , 'comment':'!'+'velocity plot options \n'},
                {'value':'1,2,1'                                          , 'comment':'!'+'bathymetry plot options \n'},
                {'value':'1'                                              , 'comment':'!'+'mesh plot option \n'},
                {'value':'0,1,481'                                        , 'comment':'!'+'mesh plot stations option \n'},
                {'value':'off'                                            , 'comment':'!'+'mesh plot element labels \n'},
                {'value':'off'                                            , 'comment':'!'+'mesh plot node lables \n'},
                {'value': plot_order                                      , 'comment':'!'+'order of nodal set for plotting straight elements \n'},
                {'value':'1'                                              , 'comment':'!'+'order of nodal set for plotting curved elements \n'},
                {'value':'1d-1,1d-1'                                      , 'comment':'!'+'zeta error tolerances: absolute, relative \n'},
                {'value':'1d-1,1d-1'                                      , 'comment':'!'+'velocity error tolerances: absolute, relative \n'},
                {'value':'1d-1,1d-1'                                      , 'comment':'!'+'barthymetry error tolerances: absolute, relative \n'},  
                {'value':'0'                                              , 'comment':'!'+'adaptive plotting option \n'},                
                {'value': zoom_box                                        , 'comment':'!'+zoom_name+' zoom box \n'},
                {'value':'7'                                              , 'comment':'!'+'figure width \n'},
                {'value':'29,33'                                           , 'comment':'!'+'timesnap range \n'},
                {'value':'/home/sbrus/Codes/dgswe/plot/work/default2.cmap', 'comment':'!'+'colormap path \n'},
                {'value':'file'                                           , 'comment':'!'+'zeta colorscale limits \n'},      
                {'value':'file'                                           , 'comment':'!'+'velocity colorscale limits \n'},    
                {'value':'auto-snap'                                      , 'comment':'!'+'bathymetry colorscale limits \n'},
                {'value':'12'                                             , 'comment':'!'+'font size \n'},    
                {'value':'cm'                                             , 'comment':'!'+'font \n'},
                {'value':'4'                                              , 'comment':'!'+'number of x ticks \n'},              
                {'value':'auto'                                           , 'comment':'!'+'number of y ticks \n'},         
                {'value':'10'                                             , 'comment':'!'+'number of c ticks \n'},      
                {'value':'2'                                              , 'comment':'!'+'number of x decimal places \n'},   
                {'value':'2'                                              , 'comment':'!'+'number of y decimal places \n'},    
                {'value':'2'                                              , 'comment':'!'+'number of c decimal places \n'},  
                {'value':'2'                                              , 'comment':'!'+'number of t decimal places \n'},      
                {'value':'png,1'                                          , 'comment':'!'+'additional file format, remove ps file option \n'},  
                {'value':'400'                                            , 'comment':'!'+'raster density \n'},   
                {'value':'0'                                              , 'comment':'!'+'movie flag \n'}        
                ] 
                
  write_file(plotting_direc+'plot.inp' ,plot_input)              

  os.chdir(plotting_direc)
  output = subprocess.Popen('./plot', shell=True, stdout=subprocess.PIPE)
  
  while output.poll() is None:
    out = output.stdout.readline()
    if out:
      print out.strip()
      
    
  subprocess.Popen('mv *.png '+save_direc, shell=True, stdout=subprocess.PIPE).communicate()[0]  
  subprocess.Popen('mv *.cscale '+save_direc, shell=True, stdout=subprocess.PIPE).communicate()[0]
  subprocess.Popen('mv plot.inp '+save_direc, shell=True, stdout=subprocess.PIPE).communicate()[0]
  subprocess.Popen('cp dgswe.inp '+save_direc, shell=True, stdout=subprocess.PIPE).communicate()[0]    
  subprocess.Popen('cp plot '+save_direc, shell=True, stdout=subprocess.PIPE).communicate()[0]    
  
