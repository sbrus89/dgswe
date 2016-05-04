function writefort14_incprec(finfile,foutfile)

[EToV,VX,B,opedat,boudat,title] = readfortfull14( strtrim(finfile) ) ;



% -- Write refinemesh to a file
%
fid = fopen(strtrim(foutfile),'w') ; 

% AGRID
fprintf(fid,'%s \n', title) ; 

% NE, NP
[NE,~] = size(EToV) ;
[NP,~] = size(VX) ;
fprintf(fid,'%d %d\n', NE, NP ) ; 

% x-y coordinate
for i = 1: NP
    fprintf(fid,'%d %25.17e %25.17e %25.17e \n', i, VX(i,1:2), B(i) ) ;  
end

% Element connectivity
for i = 1: NE
    fprintf(fid,'%d %d %d %d %d\n', i, 3, EToV(i,1:3) ) ;
end

fprintf(fid,'%d %s\n', opedat.nope, '  % Number of specified-elevation boundary segments' ) ;
fprintf(fid,'%d %s\n', opedat.neta, '  % Total number of specified-evleaction boundary nodes' ) ;
for i = 1: opedat.nope
    fprintf(fid,'%d %d %s %d \n', opedat.nvdll(i), 0, '  % open bc segment ', i ) ;
    for j = 1: opedat.nvdll(i)
        fprintf(fid, '%d \n', opdedat.nbdv(j,i) ) ; 
    end
end

fprintf(fid,'%d %s \n', boudat.nbou, '  % Number of normal-flow-boundary segments' ) ;
fprintf(fid,'%d %s \n', boudat.nvel, '  % Total number of land-boundary nodes') ;
for i = 1: boudat.nbou
   fprintf(fid,'%d %d %s %d \n', boudat.nvell(i), boudat.ibtype(i), '  % land bc segment ', i ) ;
   for j = 1: boudat.nvell(i)
       fprintf(fid, '%d \n', boudat.nbvv(j,i)) ;
   end
end

fclose(fid) ;