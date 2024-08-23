% y-difference on period:  T is fixed; f is looked for
% t = [t0 TE/2 TE]

function d = DY1(f,h,CN, t,WImet,WIop,aa)
global y000
fifi(f,h,CN);                                  warning('off','all');
[t_,y,n_] = WImet(@F_DNK,t,y000,WIop,aa);      warning('on','all');
%d = y(3,CN)-y000(CN);                                  
d = max(abs(y(3,1:2)-y000(1:2)));