function [ect_new,nelnds_new] = filter_ect(ect,nelnds,in)

ne = length(ect);
nn = length(in);

ect_tmp = zeros(ne,4);
nelnds_tmp = zeros(ne,1);
ndnums_new = zeros(nn,1);

i = 0;
for nd = 1:nn
    if (in(nd) == 1) 
        i = i +1;
        ndnums_new(nd) = i;
    end
end 

i = 0;
for el = 1:ne
   nv = nelnds(el);
   if all(in(ect(el,1:nv))) 
      i = i + 1;
      ect_tmp(i,1:nv) = ndnums_new(ect(el,1:nv));
      nelnds_tmp(i) = nv;
   end
end

ect_new = ect_tmp(1:i,:);
nelnds_new = nelnds_tmp(1:i);