function ect_new = clean_ect(ect,bnd_flag)

ne = length(ect);

ect_tmp = zeros(ne,3);

i = 0;
for el = 1:ne
   if all(bnd_flag(ect(el,:))) 
       
   else       
      i = i + 1;
      ect_tmp(i,:) = ect(el,:); 
   end
end

ect_new = ect_tmp(1:i,:);

