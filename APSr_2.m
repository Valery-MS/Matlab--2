% Almost period search,  refined appr
% y - m*n array of solution, m=numb of t-points, m=numb of coordinates
% bn - t-boundaries numbers 

function [y, T] = APSr_2(t,y0,CN,K,aa,WImet,RPmet,WIop,RPop)

tr = nan(1,K);  mif = tr;                          warning('off','all'); 
[to, y, nuf] = WImet(@RHFN,t,y0,WIop,aa,4);        warning('on','all'); 
         
[M, I] = maxk(y(:,CN),K);
Is = sort(I);
k = 1;  i = 0;
while k <= K
   S = Is(k);  
   if k<K && Is(k+1)-S==1, tL = t(S);      k = k+2;
   elseif S > 1,           tL = t(S-1);    k = k+1;
      else,                tL = t(1)-t(2); k = k+1; end
   i = i+1;
   [tr(i), mif(i)] = fminbnd(...
   @(t_) er_t(t_,RPmet,@F_DNK,t(S),y(S,:),CN,RPop,aa),...
   tL,t(S+1),optimset('TolX',eps)); end 

tr = tr(1:i);   mif = mif(1:i);
d1 = tr(2:i)  - tr(1:i-1);
d2 = d1(2:i-1)- d1(1:i-2);
j  = find( abs(d2)<2 );
if isempty(j), T  = d1(end);
else,          T  = d1(j); end