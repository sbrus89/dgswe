[ect,xy,hb,~,~,~,~] = readfort14('./fort.14_rimls');

pdeplot( xy', [], ect', 'xydata',hb, 'zdata',hb,'colormap', 'jet', 'mesh','on') ;