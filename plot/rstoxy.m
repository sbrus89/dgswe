function [x,y] = rstoxy( r, s, v )
% (r,s) --> (x,y)
%
vp = v' ; 
x = zeros(length(r),1) ;
y = zeros(length(r),1) ; 
for i = 1: length(r)
    %
    xv = -((r(i) + s(i))/2)*vp(:,1) + ...
        ((r(i) + 1)/2)*vp(:,2) + ((s(i) + 1)/2)*vp(:,3) ;
    x(i,1) = xv(1) ;
    y(i,1) = xv(2) ;
end