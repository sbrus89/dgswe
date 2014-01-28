clc
clear all

input = fopen('dgswe_backup.inp','r+');
nl = 0;
while( ~feof(input))
    nl = nl + 1;
    inp{nl,1} = fgetl(input);
end
fclose(input);

ngrids = 2;
ninp_lines = 10;

% npart(1:19) = 2:20;
% npart(20:26) = 25:5:55;

npart = [60:5:100];

% npart = [20:20 25:5:55];

n = length(npart);

for j = 3 %1:ngrids
    
    inp2 = inp;
    for line = 1:ninp_lines 
        inp2{(j-1)*ninp_lines + line,1} = inp{(j-1)*ninp_lines + line,1}(2:end);
    end
    
    input2 = fopen('dgswe.inp','w');
    for line = 1:nl
       fprintf(input2,'%s\n',inp2{line,1}); 
    end
    
    fclose(input2);

    grid = ['../grids/inlet',num2str(j),'.grd'];
   
    for i = 1:n
        
        file = fopen('partition.d','w');
        
        fprintf(file,'%s\n',num2str(npart(i)));
        fprintf(file,'%s\n',grid);
        
        !./adcprep < partition.d
        
        !./dgswe
        
    end

    
end