% Period  error depending on iv at CN coordinate of solution
% Kind of error depends on key  
            
function dy = er_ivT(t,X,T,APS,CN,aa,WImet,RPmet,WIop,RPop)
global y000   

fifi(X,CN);        % change of y000
dy = APS(t,y000,T,CN,aa,WImet,RPmet,WIop,RPop); 