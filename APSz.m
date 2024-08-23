% Almost Period Search ( zero approx, rough )
% y - m*n array of solution, m=numb of t-points, m=numb of coordinates

function [y, T] = APSz(t,iv,err,y0,CN,a,Tol,WImet,WIop) 
  
T = Inf(2,1);   
warning('off','all'); 
[to,y,nuf] = WImet(@RHFN,t,y0,WIop,a,4);        
warning('on','all');          
yi = y(:,CN);
Nopp = find(iv*yi<0,1);  % Number of 1st element opposed to y0
yi(1:Nopp) = NaN;        % Turn off value y0 from yi
  
                         % Double selection
rt = err(yi(Nopp+1:end),iv);
n_ = find(rt < Tol);             % rt numbers: rt<Tol (minimums ascending)
if isempty(n_), return;  end

rtn = rt(n_);                    % selection 1 by rt<Tol
q = 0; x = 1; cs = 0;      % cs - cч соседей      
nN = numel(n_);            % Исключаются соседние с минимумами точки -
for j = 1:nN               % тот же min, но погрешность больше
   if j < nN  &&  n_(j+1)-n_(j) == 1       %  selection 2 
      if rtn(j+1) < rtn(j), x = j+1;  end  % by condition
      cs = cs+1;                           % neighbors(n_(j+1)-n_(j)=1)
      continue,end  
  
   if x < j-cs, x = j-cs; end              % are deleting
   q = q+1; 
   T(q) = t(Nopp+n_(x));
   if q == 2, break; end
   cs = 0;  end  