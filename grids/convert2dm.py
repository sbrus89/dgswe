import pprint

filename = "/home/sbrus/data-drive/galveston/dgswe/galveston_curve.2dm"
fort14name = "/home/sbrus/data-drive/galveston/dgswe/galveston_curve.grd"
grid_name = "galveston_curve"

f = open(filename)
grid_data = f.read().split("\n")
#print grid_data

n = len(grid_data)

ect = []
xy = []
ns = []
bnd = []
bc = {}
for i in range(0,n):
  line = grid_data[i].split()
  print line
  
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
    neta = neta + len(bnd[i])
  elif bc[seg] != "0" and bc[seg] < "100":
    nbou = nbou + 1
    nvel = nvel + len(bnd[i])
    
print ne, nn
pprint.pprint(ect)
pprint.pprint(xy)
pprint.pprint(bnd)

f.close()

fort14 = open(fort14name,"w+")
fort14.write(grid_name + "\n")
fort14.write(str(ne) + " " +str(nn) + "\n")
for i in range(0,nn):
  fort14.write( "   ".join(xy[i]) + "\n")

for i in range(0,ne):
  fort14.write("   ".join(ect[i]) + "\n")
  
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
  if bc[seg] != "0" and bc[seg] < "100":
    nbndnds = len(bnd[i])
    fort14.write(str(nbndnds) + "  " +str(bc[seg]) + "\n")
    for j in range(0,nbndnds):
      fort14.write(str(bnd[i][j]) + "\n")
    
fort14.close()
