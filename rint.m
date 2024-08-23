% Root Interval
% Checking the periodicity of the functions y3,y4 is unnecessary

function AB = rint(fs,t0T,j,A0,B0,n,WImet,WIop,aa)
global y000

AB = [];
if     j == 1, fL = fs(j);   
elseif j == 2, fL = fs(j-1); 
else,          fL = fs(j-2); end

hf = (B0-fL)/10;
I = 1:2;
fifi(fL);                                            warning('off','all');
[t_,y_,n_] = WImet(@RHFN,t0T,y000,WIop,aa,4);        warning('on','all');
yL = y_(end,I)-y000(I);
f  = fL+hf;    
while f <= B0    
   fifi(f);                                           warning('off','all');
   [t_,y,n_] = WImet(@RHFN,t0T,y000,WIop,aa,4);       warning('on','all');
   yR = y(end,I)-y000(I);
   zn = yL.*yR;
   if all(zn <= 0), AB = [fL, f]; break,end
   f = f+hf;  end

if isempty(AB)
   fifi(A0);                                          warning('off','all');
   [t_,y_,n_] = WImet(@RHFN,t0T,y000,WIop,aa,4);      warning('on','all');
   yL = y_(end,I)-y000(I); 
   hf = (B0-A0)/n;
   for f = A0+hf : hf: B0
      fifi(f);                                        warning('off','all');
      [t_,y,n_] = WImet(@RHFN,t0T,y000,WIop,aa,4);    warning('on','all');
      yR = y(end,I)-y000(I);
      zn = yL.*yR;
      if all( zn <= 0), AB = [A0, f]; break,end,end,end