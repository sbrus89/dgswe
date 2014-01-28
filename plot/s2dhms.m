function [days hr min sec] = s2dhms(s)
   days = floor(s/(24*60*60));
   r = s - days*24*60*60;
   hr = floor(r/(60*60));
   r = r-hr*60*60;
   min = floor(r/60);
   r = r-min*60;
   sec = r;
      
end