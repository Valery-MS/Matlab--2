% Period  error depending on iv at CN coordinate of solution
% Kind of error depends on key  
            
function er = er_iv_1(t,X,err,CN,aa,WImet,RPmet,WIop,RPop)
global TE   y000   % T exact 

fifi(X,CN);        % change of y000
[y_, TE] = APSr(t,y000,CN,aa,WImet,RPmet,WIop,RPop); 
                                                          warning('off','all'); 
[to,y,nuf] = WImet(@RHFN,t(1)+[0 TE 2*TE],y000,WIop,aa,4); warning('on','all');       

er = err(y(3,CN),y(2,CN));                       % Abs or Rel err