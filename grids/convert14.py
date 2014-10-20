# Converts an all triangle grid so it can be read in by plot_sol_fill.m
#   Adds a 0 to the last column of the connectivity table

fort14name = "/home/sbrus/Codes/spline/nodes.out"
grdname = "/home/sbrus/Codes/spline/galveston_tri_spline.grd"

f = open(fort14name,'r')
grid_data= f.read().splitlines()
f.close()

g = open(grdname,'w')

g.write("%s\n" % grid_data[0])
g.write("%s\n" % grid_data[1])

grid_size = grid_data[1].split()

print grid_size
ne = int(grid_size[0])
nn = int(grid_size[1])

n = 2
for line in grid_data[n:nn+2]:
  g.write("%s\n" % line)
  n = n+1
  
for line in grid_data[n:n+ne]:
  line = line + "       0     0      0"
  g.write("%s\n" % line )
  n = n+1
  
for line in grid_data[n:]:
  g.write("%s\n" % line )
  
  
g.truncate()
g.close

