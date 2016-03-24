clc
clear all
close all


p =2;

knots = [0 0 0 1 2 3 4 4 5 5 5];
% knots = [0 1 2 3 4 5 6 7];

nk = length(knots);
n = nk-(p+1);

xi = linspace(knots(1),knots(end),1000)';

N = zeros(length(xi),n+1,p+1);

for i = 1:nk-1
    N(:,i,1) = (xi >= knots(i) & xi < knots(i+1));
end


for j = 2:p+1
    pj = j-1;
    for i = 1:n
        
        if (knots(i+pj)-knots(i)) == 0
            a = 0;
        else            
            a = (xi-knots(i))/(knots(i+pj)-knots(i));
        end
        
        if (knots(i+pj+1)-knots(i+1)) == 0 
            b = 0;
        else            
            b = (knots(i+pj+1)-xi)/(knots(i+pj+1)-knots(i+1));
        end
        
        N(:,i,j) = a.*N(:,i,j-1) + b.*N(:,i+1,j-1);
        
    end
end

xy = [0    0;
      0    1;
      1    -.75;
      1    .5;
      2    .5;
      2.5  -.75;
      3    .25;
      3.5  .25];
  
cx = xi*0;
cy = xi*0;
for i = 1:n
  cx = cx + N(:,i,3)*xy(i,1);
  cy = cy + N(:,i,3)*xy(i,2);
end


figure
hold all
for i = 1:n
    plot(xi,N(:,i,3))
end
axis image

figure
hold on
plot(xy(:,1),xy(:,2),'ro')
plot(xy(:,1),xy(:,2),'k')
plot(cx(1:end-1),cy(1:end-1))
