% Almost period search,  refined appr
% y - m*n array of solution, m=numb of t-points, m=numb of coordinates
% bn - t-boundaries numbers

function [y, T] = APSr_1(t,X,err,y0,CN,dt,NoJ,dn,a,WImet,RPmet,WIop,RPop)
global J bn ert Ts_ TG TT    % грубое точное

if isempty(J),   J = 0;   end
J = J+1;
warning('off','all'); 
[to, y, nuf] = WImet(@RHFN,t,y0,WIop,a,4);        
warning('on','all');          

if numel(bn) == 1
   rt = err(y(:,CN),X);   rt(1) = Inf;
   [m, N] = min( rt ); 
   
   if bn == 1,  T = t(N);
   else
      rt(N-1:N+1) = Inf;
      [m, m2] = min(rt);   MM = [N m2];
      TG(:,J) = t(MM);
      k = 0; 
      for M = MM
         tM = t(M); 
         if y(M-1) > y(M+1), t1 = tM-dt;  t2 = tM;
         else,               t1 = tM;     t2 = tM+dt;   end
         k = k+1;          
         [TT(k,J), ert(k,J)] = fminbnd(...
         @(t) er_t(t,X,err,RPmet,@F_DNK,tM,y(M,:),CN,RPop,a),...
         t1,t2,optimset('TolX',eps)); end 
      [mi, K] = min(ert(:,J));
      T = TT(K,J); 
      N = MM(K); end 
                    
   Ts_(J) = T;
   if J > NoJ && all(abs(Ts_(J-NoJ+1:J)-Ts_(J-NoJ:J-1)) < 1 )
      bn = N-dn : min(N+dn,numel(t)); end
else    
   rt = err(y(bn,CN),X);
   [m, N] = min(rt);             
   M = N+bn(1)-1;   tM = t(M);

   if y(M-1) > y(M+1), t1 = tM-dt;  t2 = tM;
   else,               t1 = tM;     t2 = tM+dt;   end
                 
   [T, ert(1,J)] = fminbnd(...
                   @(t) er_t(t,X,err,RPmet,@F_DNK,tM,y(M,:),CN,RPop,a),...
                   t1,t2,optimset('TolX',eps)); 
   Ts_(J) = T; end