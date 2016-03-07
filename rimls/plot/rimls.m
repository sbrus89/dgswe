clear all
close all
clc

path(path,'/home/sbrus/Codes/SomeDGMaterials')
path(path,'/home/sbrus/Codes/SomeDGMaterials/devscript')
path(path,'/home/sbrus/Codes/SomeDGMaterials/ModNodeADCIRCDG')
path(path,'/home/sbrus/Codes/SomeDGMaterials/NodalDG')

direc = '/home/sbrus/Codes/dgswe/grids/';
name = 'inlet1.grd';

% direc = '/home/sbrus/Codes/rimls/';
% name = 'beaufort.grd';

direc2 = '/home/sbrus/Codes/rimls/';

[ect,xy,hb,~,~,~,~] = readfort14([direc,name]);

nn = length(xy);
ne = length(ect);

nepn = zeros(1,nn);
for el = 1:ne
    n1 = ect(el,1);
    n2 = ect(el,2);
    n3 = ect(el,3);
    
    nepn(n1) = nepn(n1) + 1;   
    nepn(n2) = nepn(n2) + 1;    
    nepn(n3) = nepn(n3) + 1;
end
mnepn = max(nepn);

nepn = zeros(1,nn);
epn = zeros(mnepn,nn);
for el = 1:ne
    n1 = ect(el,1);
    n2 = ect(el,2);
    n3 = ect(el,3);
    
    nepn(n1) = nepn(n1) + 1;
    epn(nepn(n1),n1) = el;
    
    nepn(n2) = nepn(n2) + 1;
    epn(nepn(n2),n2) = el;
    
    nepn(n3) = nepn(n3) + 1;
    epn(nepn(n3),n3) = el;
    
end


[h,hmax,hmin] = grid_size(ne,ect,xy);

fid = fopen([direc2,'edge_nodes.d']);
N = fscanf(fid,'%g %g',2) ;
nends = N(1)*N(2);
ends = fscanf(fid,'%g %g %g \n', [3 nends])' ;

fid4 = fopen([direc2,'interior_nodes.d']);
N = fscanf(fid4,'%g %g',2) ;
ninds = N(1)*N(2);
inds = fscanf(fid4,'%g %g %g \n', [3 ninds])' ;

fid2 = fopen([direc2,'centers.d']);
ne = fscanf(fid2,'%g',1) ;
cnds = fscanf(fid2,'%g %g %g \n', [3 ne])' ;

fid3 = fopen([direc2,'normals.d']);
ne = fscanf(fid3,'%g',1) ;
n = fscanf(fid3,'%g %g %g \n', [3 ne])' ;


kdtreeobj = KDTreeSearcher(xy);
ckdtreeobj = KDTreeSearcher(cnds);

xpts = vertcat(xy(:,1),ends(:,1),inds(:,1));
ypts = vertcat(xy(:,2),ends(:,2),inds(:,2));
hpts = vertcat(hb,ends(:,3),inds(:,3));

figure(1)
ect2 = delaunay(xpts,ypts);
nodes = [xpts ypts];
pdeplot( nodes', [], ect2', 'xydata',hpts,'zdata',hpts, 'colormap', 'jet', 'mesh','off') ;

nnr = length(xpts);

max_iter = 10;
threshold = 1e-4;
sigma_r = 0.5;
sigma_n = 5; % 0.5-1.5
range = 4; % 1.5-4


hbv = zeros(nnr,1);
x = zeros(1,3);
for ipt = 1:nnr
    disp([num2str(ipt),'/',num2str(nnr)])    
    
    x = [xpts(ipt) ypts(ipt) hpts(ipt)];
    
    grad_f = ones(1,3);
    f = 1;
    
    nd = knnsearch(kdtreeobj,x(1:2));
    
    lh = range*h(epn(1,nd));    

    idx = rangesearch(ckdtreeobj,x,lh);
    
    nneigh = length(idx{1});

    while norm(f*grad_f) > threshold
        it = 0;
        
        alpha = ones(nneigh,1);
        alpha_old = zeros(nneigh,1);
        
        while it < max_iter
%         while max(abs(alpha/sum(alpha) - alpha_old/sum(alpha_old))) > threshold %|| it < max_iter            
            
            sumW = 0;
            sumGw= zeros(1,3);
            sumF = 0;
            sumGF = zeros(1,3);
            sumN = zeros(1,3);
            
            for cpt = 1:nneigh
                p = cnds(idx{1}(cpt),:);
                np = n(idx{1}(cpt),:);                               
                px = x - p;
                
                fx = dot(px,np);
                
                if it > 0 
                    alpha(cpt) = exp(-((f-fx)/(sigma_r*lh))^2)*exp(-(norm(np-grad_f)/sigma_n)^2);
                else
                    alpha(cpt) = 1;
                end
                
                w = alpha(cpt)*(1-norm(px)^2/lh^2)^4;
                
                grad_w = alpha(cpt)*grad_phi(px,lh);

                sumW = sumW + w;
                sumGw = sumGw + grad_w;
                sumF = sumF + w*fx;
                sumGF = sumGF + grad_w*fx;
                sumN = sumN + w*np;
            end
            
            f = sumF/sumW;
            grad_f = (sumGF - f*sumGw + sumN)/sumW;
       
            it = it +1;

            alpha_old = alpha;            
            
        end
       
        x = x - f*grad_f;
       
    end
    
    xpts(ipt) = x(1);
    ypts(ipt) = x(2);
    hpts(ipt) = x(3);
    
end

figure(2)
ect2 = delaunay(xpts,ypts);
nodes = [xpts ypts];
pdeplot( nodes', [], ect2', 'xydata',hpts,'zdata',hpts, 'colormap', 'jet', 'mesh','off') ;
