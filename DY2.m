% comparison by half-period 
% CN - only 1 or 2; CNG = only 3 or 4
% t = [t0 TE/2]

function d = DY2(f,h,CN, t,WImet,WIop,aa)
global y000 

fifi(f,h,CN);                                           warning('off','all');
[t_,y,n_] = WImet(@F_DNK,t,y000,WIop,aa);               warning('on','all');
% d = y(end,3);        
d = max(abs(y(end,3:4)));  