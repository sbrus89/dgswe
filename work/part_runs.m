clc
clear all

input = fopen('dgswe_backup.inp','r+');
nl = 0;
while( ~feof(input))
    nl = nl + 1;
    inp{nl,1} = fgetl(input);
end
fclose(input);

ngrids = 1;
ninp_lines = 12;

npart = 1:2;

n = length(npart);

for j = 1:ngrids
    
    for i = 1:n
        
        inp2 = inp;
        for line = 1:ninp_lines -2
            inp2{(j-1)*ninp_lines + line,1} = inp{(j-1)*ninp_lines + line,1}(2:end);
        end
        line = line + 1;
        inp2{(j-1)*ninp_lines + line,1} = num2str(npart(i));
        
        input2 = fopen('dgswe.inp','w');
        for line = 1:nl
            fprintf(input2,'%s\n',inp2{line,1});
        end
        
        fclose(input2);
        
        !./dgswe
        
    end
    
    
end