function [sol,t,nsnap_max] = read_solution(direc,fname,nsnap_read)

fid = fopen([direc,fname]);

header = 1;
while header == 1
    line = fgetl(fid);
    disp(line)
    if length(line) > 10 && strcmp(line(1:10),'!!!!!!!!!!');
        header = 0;
    end
end

line = fscanf(fid,'%d',[3,1])';
disp(line)
ndof = line(1);
ne = line(2);
nsnap = line(3);

if nsnap_read > nsnap   
   nsnap_read = nsnap; 
end

sol = zeros(ndof,ne,nsnap_read);
t = zeros(nsnap_read,1);

snap = 0;
while ~feof(fid) && snap < nsnap_read
    snap = snap + 1;
    t(snap) = fscanf(fid,' %g ', 1);     
    sol(:,:,snap) = fscanf(fid,' %g ', [ne ndof])'; 
end

if snap < nsnap_read
   disp('warning: fewer time snaps read in than expected') 
end

nsnap_max = snap;
end

