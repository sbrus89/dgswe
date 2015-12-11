clc
close all
clear all

fid = fopen('../output/bathy.d','r') ;

i = 1;
while ~feof(fid)
    Val(i,1:3) = fscanf(fid,'%g %g %g \n', 3 ) ;
    i = i + 1;
end

xy(:,1:2) = Val(:,1:2);
hb = Val(:,3);

ect = delaunay(xy(:,1),xy(:,2));
pdeplot( xy', [], ect', 'xydata',hb, 'zdata',hb, 'colormap', 'jet','mesh','on') ;


figure

hb2 = 10 - 5*cos(2*pi/500*xy(:,2));
pdeplot( xy', [], ect', 'xydata',hb2,'zdata',hb2, 'colormap', 'jet','mesh','on') ;