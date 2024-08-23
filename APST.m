% Almost period search:  for a given T, find y0
% y - m*n array of solution, m=numb of t-points, m=numb of coordinates

function dy = APST(t0,f,T,aa,WImet,WIop)
global y000 Y

CN = 1;
fifi(f,CN);                                             warning('off','all'); 
[to,Y,nuf] = WImet(@RHFN,t0+[0 T/2 T],y000,WIop,aa,4);  warning('on','all'); 
         
dy = max(abs(Y(3,1:1)-y000(1:1)));