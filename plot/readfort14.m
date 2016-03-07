function [EToV,VX,B,nelnds,opedat,boudat,title] = readfort14( finame ) 
% Read ADCIRC grid file fort.14

tic
if ( nargin == 0 ) 
  finputname = 'fort.14';
else
  finputname = finame ;   
end

fid = fopen(finputname) ;

agrid = fgetl(fid) ;
disp(agrid) ;
title = agrid ; 

N = fscanf(fid,'%g %g',2) ;

% curved_grid = 0;
% if N(2)/N(1) > 3 % might need to be 2 for curved tri grids
%    curved_grid = 1 
% end

Val = fscanf(fid,'%d %g %g %g \n', [4 N(2)])' ;

iv = sort(Val(:,1)) ;
Val = Val(iv,:) ; 

% if curved_grid == 1
%     idx = fscanf(fid,'%d %d %d %d %d %d %d %d %d %d %d \n', [11 N(1)])' ;
% else
% %     idx = fscanf(fid,'%d %d %d %d %d %d \n', [6 N(1)])' ;
%     idx = fscanf(fid,'%d %d %d %d %d \n', [5 N(1)])' ;
% end

idx = zeros(N(1),11);
for i = 1:N(1)
    line = fgetl(fid);
    v = sscanf(line,'%d',11)';
    idx(i,1) = v(1);
    idx(i,2) = v(2);
    idx(i,3:v(2)+2) = v(3:v(2)+2);
end

iv = sort(idx(:,1)) ;
idx = idx(iv,:) ; 

% Arrange it to a Nodal DG input
VX = Val(:,2:3) ;
B  = Val(:,4) ;
% if curved_grid == 1
%     EToV = idx(:,3:11) ;
% else
% %     EToV = idx(:,3:6) ;
%     EToV = idx(:,3:5);
% end
EToV = idx(:,3:end);
nelnds = idx(:,2);

% Read in boundary
% Open boundary
msgline = fgetl(fid) ;
nope = sscanf(msgline,'%d %*s') ;

msgline = fgetl(fid) ; 
neta = sscanf(msgline,'%d %*s') ;

nvdll = zeros(nope,1) ;
ibtypee = zeros(nope,1) ;
nbdv = zeros(neta,nope) ; 
% tic
for i = 1: nope
   msgline = fgetl(fid) ;
   
   [varg] = sscanf(msgline,'%d %*s \n') ;
   nvdll(i) = varg ;
   ibtypee(i) = 0 ;
   
   %
   % for k = 1: nvdll(i)
   %    
   %    % % nbdv(k,i) = fscanf(fid,'%d \n')  
   %    % % fscanf(fid,'%d \n')
   %    msgline = fgetl(fid) ; 
   %    
   %    nbdv(k,i) = str2num(msgline) ; 
   % end
   %
   % Nov 25, 2012, improve reading efficiency
   nbdv(1:nvdll(i),i) = fscanf(fid,'%d \n', nvdll(i) ) ;  
end
% toc

opedat.nope = nope ; 
opedat.neta = neta ;
opedat.nvdll = nvdll ;
opedat.ibtypee = ibtypee ;
opedat.nbdv = nbdv ; 

% land boundary
msgline = fgetl(fid) ;
nbou = sscanf(msgline,'%d %*s') ;

msgline = fgetl(fid) ; 
nvel = sscanf(msgline,'%d %*s') ;

nvell = zeros(nbou,1) ;
ibtype = zeros(nbou,1) ;
nbvv = zeros(nvel,nbou) ;
% tic
for i = 1: nbou
    msgline = fgetl(fid) ;
    
    [varg,nag] = sscanf(msgline,'%d %d %*s \n') ;
    nvell(i) = varg(1) ;
    ibtype(i) = varg(2) ;
    
    switch ( ibtype(i) )
        case {0,1,2,10,11,12,20,21,22,30,60,61,101,102}
            %for k = 1: nvell(i)
            %   msgline  = fgetl(fid) ;
            %   nbvv(k,i) = str2num(msgline) ;
            %end
            % Nov 15, 2012, improve reading efficiency
            nbvv(1:nvell(i),i) = fscanf(fid,'%d \n', nvell(i) ) ;
        otherwise
            msgline = fgetl(fid) ;
    end
end
% toc

% land boundary
boudat.nbou = nbou ;
boudat.nvel = nvel ;
boudat.nvell = nvell ;
boudat.ibtype = ibtype ;
boudat.nbvv = nbvv ;  


fclose(fid) ;

return

