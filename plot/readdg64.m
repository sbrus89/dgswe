function [Q,tt,nsnap,Nel,ntot] = readdg64( nsnapi,dir )
% Read ADCIRC DG.64 (uH, vH) 

% read output file
% tic
fid1 = fopen([dir,'DG.64'],'r') ;

% Discard
fline = fgetl(fid1) ;
disp('Read DG.64:') ;
disp(fline) ;

v2 = fscanf(fid1,'%f  %d %f %d %d',5) ;

% get the number of elements
fid2 = fopen([dir,'fort.14']) ;
agrid = fgetl(fid2) ;
N = fscanf(fid2,'%g %g',2) ;
fclose(fid2) ;

Nel = N(1) ;

% get the number of total nodes
ntot = Nel*v2(2) ;

if ( nargin == 0 ) 
   nsnapi = v2(1) ; 
end
% Snapshots to be recorded
nsnap = 0 ;
tt = zeros(nsnapi,2) ;
Q  = zeros(ntot, 2, nsnapi) ;
for i = 1: nsnapi
    % if ( feof(fid1) )
    %    break ;
    % end
    % tt(i,:) = fscanf(fid1,'%f  %f',2) ;
    
    tval =  fscanf(fid1,'%f  %f',2) ;
    fed = 0 ;
    while ( isempty(tval) )
        fed = feof(fid1) ;
        if ( fed )
            break ; 
        else
            tval =  fscanf(fid1,'%f  %f',2) ;
        end
    end
    if ( ~fed )
        tt(i,:) = tval ; 
    else
        break ;
    end
    
    %
    % for n = 1: ntot
    %    val = fscanf(fid1,'%d %f %f',3) ;
    %    Q(n,1:2,i) = val(2:3) ;
    % end
    %
    
    % Nov 15, 2012, improve efficiency
    val = fscanf(fid1, '%d %f %f \n', [3, ntot])' ; 
    Q(:,1:2,i) = val(:,2:3) ; 
    
    nsnap = nsnap + 1 ;
    fprintf('%c','.') ; 
    if ( mod(i,80) == 0 )
        fprintf('\n') ; 
    end
end
fprintf('\n') ;

fclose(fid1) ;
% toc