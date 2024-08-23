% Almost period search,  refining of period 
% y - m*n array of solution, m=numb of t-points, n=numb of coordinates
% T  is not used here  

function [dy TE] = APSr(t,f,h,CN,aa,WImet,RPmet,WIop,RPop)
global y000 TE Yma                % Yma  - from er_t

nt = numel(t);  
K  = 6;     tr = nan(K,1);  Ys = nan(K,4);   
fifi(f,h,CN);                                      warning('off','all'); 
[to, y, nuf] = WImet(@F_DNK,t,y000,WIop,aa);       warning('on','all'); 
         
[M,I] = maxk(y(:,CN),K); % search for К maxima from the coordinates of y(:,CN)
k = 1;  i = 0;
while k <= K
   J = I(k);  
   if     J == 1,                 tL = t(J);   tR = t(J+1);   k = k+1;
   elseif J == nt,                tL = t(J-1); tR = t(J);     k = k+1;
   elseif k<K && abs(I(k+1)-J)<2, tL = t(J);   tR = t(J+1);   k = k+2;
   else,                          tL = t(J-1); tR = t(J+1);   k = k+1; end
   i = i+1;
   tr(i) = fminbnd(@(t_) er_t(t_,RPmet,@F_DNK,t(J),y(J,:),CN,RPop,aa),...
           tL,tR,optimset('TolX',eps)); 
           % - точка локального максимума ф-ции y(t,CN) на [tL,tR]
   Ys(i,:) = Yma; end   % - локальный максимум этой ф-ции

[MM,II] = maxk(Ys(:,CN),2);
dy = MM(1)-MM(2);
TE = abs(tr(II(1))-tr(II(2)));
