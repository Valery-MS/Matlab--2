% Calc fi from energy integral
% CN - only 1 or 2

function fifi(f,h,CN)
global a b y000

co = cos(f);  
A  = sin(f);        B = b-co;   
C  = a*(1-co)+b-h;  S = sqrt(A*A+B*B);
if abs(C) <= S         %, waitfor(errordlg('fifi: |C| > S'));end
   y000(CN)   = f; 
   y000(3-CN) = asin(C/S) - atan(B/A);
end

%{
if abs(C) > S 
  e=1e-10;
  x=1-co;  a2=a*a;
  hmi=a+b-(a*b+a/b+b/a)/2;
  xc=a*(h-b)+b; D=sqrt(2*a*b*(h-hmi)); 
  xmi=(xc-D)/a2; xma=(xc+D)/a2;
  if x<=xmi || x>=xma
     [xmi x xma],end
  c=a*x+b-h; s=sqrt((b-1)^2+2*b*x);
  if abs(C-c)>r || abs(S-s)>e
     [C c S s],end
end
%}