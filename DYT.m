% Difference on period:  f is fixed; T is looked
% t = [t0 TE/2 TE]
             
function d = DYT(T,CN,WImet,WIop,aa)
global y000
                                                       warning('off','all');
[t_,y,n_] = WImet(@F_DNK,[0 T/2 T],y000,WIop,aa);      warning('on' ,'all');
%d = y(3,CN)-y000(CN);                                  
d = abs(y(3,CN)-y000(CN));