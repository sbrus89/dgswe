load ../fort.64

close all

figure
pdeplot( VX', [], EToV','xydata',HB,'zdata',HB,'xystyle','interp','colormap','jet','mesh','on') ;
xlabel('x')
ylabel('y')
title('hb') 
 
figure
pdeplot( VX', [], EToV','xydata',fort(:,1),'zdata',fort(:,1),'xystyle','interp','colormap','jet','mesh','on') ;
xlabel('x')
ylabel('y')
title('ddx')

figure
pdeplot( VX', [], EToV','xydata',fort(:,2),'zdata',fort(:,2),'xystyle','interp','colormap','jet','mesh','on') ;
xlabel('x')
ylabel('y')
title('ddy')