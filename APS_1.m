% Almost period search
% y - m*n array of solution, m=numb of t-points, m=numb of coordinates
% T - [T1, T2, T2-T1]
% nmi - period number for which the minimum relative error was obtained

function [y, T, nmi] = APS_1(t,iv,err,y0,CN,h,nt,a,Tol,kT,...
                          WImet,RPmet,SPmet,WIop,RPop,SPop)
global yTA   % calcs in solt с.164:  Tj = fminbnd(@(t) -metR(t,... )

y0(CN) = iv;
dt = t(2)-t(1);   dt_ = dt*0.7; 
T  = Inf(kT,1);
ers = T;                          % Table of Relative Tolerances
warning('off','all'); 
[to,y,nuf] = WImet(@RHFN,t,y0,WIop,a,4);        
warning('on','all');          
yi = y(:,CN);
Nopp = find(iv*yi<0,1);  % Number of 1st element opposed to y0
yi(1:Nopp) = NaN;        % Turn off value y0 from yi
      
                         % Double selection
rt = err(yi(Nopp+1:end),iv);
n_ = find(rt < Tol);             % rt numbers: rt<Tol (minimums ascending)
rtn = rt(n_);                    % selection 1 by rt<Tol

q = 0; x = 1; cs = 0;  I = n_;   % cs - cч соседей      
nN = numel(n_);            % Исключаются соседние с минимумами точки -
for j = 1:nN               % тот же min, но погрешность больше
   if j < nN  &&  n_(j+1)-n_(j) == 1       %  selection 2 
      if rtn(j+1) < rtn(j), x = j+1;  end  % by condition
         cs = cs+1;                        % neighbors(n_(j+1)-n_(j)=1)
      continue,end              
   if x < j-cs, x = j-cs; end              % are deleting
   q = q+1;  cs = 0;
   I(q) = n_(x);  end
Ms = Nopp+I(1:q);     
           
j = 0; r = 0; 
ert = nan(q,1);
while r < kT  &&  j < q    
   j = j+1;  M = Ms(j);   tM = t(M);
   if r > 0                  % search for a second (or more) AP
      dr = tM./T(1:r);
      if any(abs(dr-round(dr))<0.05), continue,end,end
        
   if     M == nt,           t1 = tM-dt_;  t2 = tM+dt; 
   elseif yi(M-1) > yi(M+1), t1 = tM-dt_;  t2 = tM;
   else                      t1 = tM;      t2 = tM+dt_;   end
                 
   [Tj ert(j)]= fminbnd(@(t)er_t(t,iv,err,RPmet, @F_DNK,tM,y(M,:),CN,RPop,a),...
                t1,t2,optimset('TolX',eps));   
 
   if ert(j) < Tol, r=r+1;  T(r) = Tj;  ers(r) = ert(j); end
end
      
%if
if ers(1)<ers(2), nmi = 1; else nmi = 2; end 