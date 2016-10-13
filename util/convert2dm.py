import pprint

#filename = "/home/sbrus/data-drive/galveston/dgswe/quad2_spline/galveston2.2dm"
#fort14name = "/home/sbrus/data-drive/galveston/dgswe/quad2_spline/galveston2_spline_plot.grd"
#grid_name = "galveston2_spline_plot"

#filename = "/home/sbrus/Codes/spline/galveston2_curve.2dm"
##fort14name = "/home/sbrus/Codes/spline/galveston2_spline.grd"
#fort14name = "/home/sbrus/data-drive/galveston/dgswe/quad2_spline_channel/galveston2_spline.grd"
#grid_name = "galveston2_spline"

#filename = "/home/sbrus/Codes/spline/galveston_tri_curve.2dm"
#fort14name = "/home/sbrus/Codes/spline/galveston_tri_curve.grd"
#grid_name = "galveston_tri_spline"

#filename = "/home/sbrus/Codes/error/converge_quad3.2dm"
#fort14name = "/home/sbrus/Codes/error/converge_quad3.grd"
#grid_name = "converge_quad3"

#filename = "/home/sbrus/data-drive/converge_quad/converge_quad4.2dm"
#fort14name = "/home/sbrus/data-drive/converge_quad/converge_quad4.grd"
#grid_name = "converge_quad4"

filename = "/home/sbrus/data-drive/galveston_spline_flux_fix/grids/galveston_straight_quad_oob.2dm"
fort14name = "/home/sbrus/data-drive/galveston_spline_flux_fix/grids/galveston_quad.grd"
grid_name = "galveston_quad_flux"

f = open(filename)
grid_data = f.read().splitlines()
#print grid_data

n = len(grid_data)

curved = 0

ect = []
xy = []
ns = []
bnd = []
bc = {}
for i in range(0,n):
  line = grid_data[i].split()
  #print line
  
  data = []
  
  if not line:  # check if empty 
    print "line is empty"
  else:
    
    if line[0] == "E4Q":
      data.append(line[1])
      data.append("4")
      data.append(line[2])
      data.append(line[3])
      data.append(line[4])
      data.append(line[5])    
      ect.append(data)
    elif line[0] == "E3T":
      data.append(line[1])
      data.append("3")
      data.append(line[2])
      data.append(line[3])
      data.append(line[4])
      data.append("0")
      ect.append(data)
    elif line[0] == "ND":
      data.append(line[1])
      data.append(line[2])
      data.append(line[3])
      data.append(line[4])
      xy.append(data)
    elif line[0] == "E9Q":
      data.append(line[1])
      data.append("9")
      data.append(line[2])
      data.append(line[3])
      data.append(line[4])
      data.append(line[5])  
      data.append(line[6])
      data.append(line[7])
      data.append(line[8])
      data.append(line[9])
      data.append(line[10])      
      ect.append(data)    
      
      curved = 1
    elif line[0] == "E6T":
      data.append(line[1])
      data.append("6")
      data.append(line[2])
      data.append(line[3])
      data.append(line[4])
      data.append(line[5]) 
      data.append(line[6]) 
      data.append(line[7]) 
      data.append("0")
      data.append("0")
      data.append("0")      
      ect.append(data)  
      
      curved = 1
    elif line[0] == "NS":
      if int(line[-2]) > 0:
        for j in range(1,len(line)):
	  ns.append(line[j])
      else: 
        for j in range(1,len(line)-1):
          ns.append(str(abs(int(line[j]))))
        bnd.append(ns)
        ns = []
    elif line[0] == "BC_VAL":
      bc[line[2]] = line[3]
       
       
ne = len(ect)
nn = len(xy)
nbnd = len(bnd)
tbndnds = 0
nope = 0
neta = 0
nbou = 0
nvel = 0
for i in range(0,nbnd):
  seg = str(i+1)
  tbndnds = tbndnds + len(bnd[i])
  if bc[seg] == "0":
    nope = nope + 1
    if curved: 
      neta = neta + (len(bnd[i])+1)/2
    else:
      neta = neta + len(bnd[i])
  elif bc[seg] != "0" and int(bc[seg]) < 200:
    nbou = nbou + 1
    if curved:
      nvel = nvel + (len(bnd[i])+1)/2
    else:
      nvel = nvel + len(bnd[i])
    
print ne, nn
#pprint.pprint(ect)
#pprint.pprint(xy)
#pprint.pprint(bnd)
print nope,nbou,nbnd

f.close()

fort14 = open(fort14name,"w+")
fort14.write(grid_name + "\n")
fort14.write(str(ne) + " " +str(nn) + "\n")
for i in range(0,nn):
  fort14.write( "   ".join(xy[i]) + "\n")

for i in range(0,ne):
  fort14.write("   ".join(ect[i]) + "\n")
  
  
if curved: # need to do every other node 
  fort14.write(str(nope)       + " = Number of open boundaries\n")  
  fort14.write(str(neta) + " = Total number of open boundary nodes\n")
  for i in range(0,nbnd):
    seg = str(i+1)
    if bc[seg] == "0":
      nbndnds = (len(bnd[i])+1)/2    
      fort14.write(str(nbndnds) + "  " +str(bc[seg]) + "\n")
      for j in range(0,nbndnds):
        fort14.write(str(bnd[i][2*j]) + "\n")      
        
        
  fort14.write(str(nbou) + " = Number of land boundaries\n")
  fort14.write(str(nvel) + " = Total number of land boundary nodes\n")
  for i in range(0,nbnd):
    seg = str(i+1)
    if bc[seg] != "0" and int(bc[seg]) < 200:
      nbndnds = (len(bnd[i])+1)/2
      fort14.write(str(nbndnds) + "  " +str(bc[seg]) + "\n")
      for j in range(0,nbndnds):
        fort14.write(str(bnd[i][2*j]) + "\n")

else:
  fort14.write(str(nope) + " = Number of open boundaries\n")  
  fort14.write(str(neta) + " = Total number of open boundary nodes\n")
  for i in range(0,nbnd):
    seg = str(i+1)
    if bc[seg] == "0":
      nbndnds = len(bnd[i])    
      fort14.write(str(nbndnds) + "  " +str(bc[seg]) + "\n")
      for j in range(0,nbndnds):
        fort14.write(str(bnd[i][j]) + "\n")      
        
        
  fort14.write(str(nbou) + " = Number of land boundaries\n")
  fort14.write(str(nvel) + " = Total number of land boundary nodes\n")
  for i in range(0,nbnd):
    seg = str(i+1)
    if bc[seg] != "0" and int(bc[seg]) < 200:
      nbndnds = len(bnd[i])
      fort14.write(str(nbndnds) + "  " +str(bc[seg]) + "\n")
      for j in range(0,nbndnds):
        fort14.write(str(bnd[i][j]) + "\n")
    
fort14.close()
