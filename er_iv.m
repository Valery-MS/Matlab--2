% Period  error depending on iv at CN coordinate of solution
% Kind of error depends on key
            
function er = er_iv(t,X,err,CN,dt,NoJ,dn,K,a,WImet,RPmet,WIop,RPop)
global TE   y000   % T exact 

fifi(X,CN);        % change of y000
[y_, TE] = APSr(t,X,err,y000,CN,dt,NoJ,dn,K,a,WImet,RPmet,WIop,RPop); 
                                                        warning('off','all'); 
[to,y,nuf] = WImet(@RHFN,[t(1) t(1)+TE],y000,WIop,a,4); warning('on','all');       

er = err(y(end,CN),X);                       % Abs or Rel err