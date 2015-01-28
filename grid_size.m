function [h hmax hmin] = grid_size(ne,ect,xy)

    h = zeros(ne,1);
    
    for el = 1:ne
       x1 = xy(ect(el,1),1);
       x2 = xy(ect(el,2),1);
       x3 = xy(ect(el,3),1);
       
       y1 = xy(ect(el,1),2);
       y2 = xy(ect(el,2),2);
       y3 = xy(ect(el,3),2);
       
       a = sqrt((x1-x2)^2 + (y1-y2)^2);
       b = sqrt((x2-x3)^2 + (y2-y3)^2);
       c = sqrt((x3-x1)^2 + (y3-y1)^2);
       
       s = .5*(a+b+c);
       
       r = sqrt((s-a)*(s-b)*(s-c)/s);
       
       h(el) = 2*r;
       
    end
    
    hmax = max(h);
    hmin = min(h);
    
end