function [ze,tt,nsnap,Nel,ntot] = readdg63( nsnapi,dir )
% Read ADCIRC DG.63 (water level from geoid)

tic
% read output file
fid1 = fopen([dir,'DG.63'],'r') ;

% Discard
fline = fgetl(fid1) ;
disp('Read in DG.63:') ;
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

% Snapshot to be recorded
nsnap = 0 ;
zedum = zeros(ntot,2) ;

tt = zeros(nsnapi,2) ;
ze = zeros(ntot,2,nsnapi)  ;

for i = 1: nsnapi
    %
    %
    tval = fscanf(fid1,'%f  %f \n',2) ;
    fed = 0 ;
    while ( isempty(tval) )
        fed = feof(fid1) ;
        if ( fed )
            break ; 
        else
          tval = fscanf(fid1,'%f  %f \n',2) ;
        end
    end
    
    if ( ~fed )
        tt(i,:) = tval ; 
    else
        break ; 
    end
    
    %
    % for n = 1: ntot
    %    msgline = fgetl(fid1) ;
    %    
    %    valn = str2num(msgline) ;
    %    [~,ic] = size(valn) ;
    %    
    %    zedum(n,1:ic-1) = valn(2:ic) ;
    %    zedum(n,ic:2) = -99999 ; 
    % end    
    % ze(1:ntot,1:2,i) = zedum ; 
    %
    
    % For Tim Stitt
    % Nov 15, 2012, improv efficiency
    %-----------------------------------------
    % msgline = fgetl(fid1) ;
    % valn = str2num(msgline) ;  %#ok<ST2NM>
    % [~,ic] = size(valn) ;
    %
    % zedum(1,1:ic-1) = valn(2:ic) ;
    % zedum(1,ic:2) = -99999 ;
    %
    % fmtread = '%d %g %g \n' ;
    % if ( ic == 2 ) 
    %    fmtread = '%d %g \n' ;
    % end
    %
    % val = fscanf(fid1, fmtread, [ic ntot-1]) ;
    % zedum(2:ntot,1:ic-1) = val(2:ic,:)' ;
    % zedum(2:ntot,ic:2) = -99999 ;
    %
    % ze(:,1:2,i) = zedum ; 
    %------------------------------------------
    
    % 
    % Nov 15, 2012, improv efficiency
    %-----------------------------------------  
    fmtread = '%d %g %g \n' ;
    
    val = fscanf(fid1, fmtread, [3 ntot]) ;
    zedum(1:ntot,1:2) = val(2:3,:)' ;
    
    ze(1:ntot,1:2,i) = zedum ; 
    %------------------------------------------
    
    nsnap = nsnap + 1 ;
    
    fprintf('%c','.') ; 
    if ( mod(i,80) == 0 ) 
        fprintf('\n') ;
    end
end
fprintf('\n') ;

fclose(fid1) ;
toc