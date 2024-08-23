% Period error depending on t at CN coordinate of solution
% CN - Chaning Initial Value Index
% (tM,yM) - extremum point in discret solution y(:,CN)
% Yma - y in local max point   

function ymi = er_t(ti,met,F,tM,yM,CN,op,aa)
global Yma

if ti == tM, Yma = yM;
else      
   op.InitialStep = (ti-tM)*0.5;
   warning('off','all');
   [t_, y] = met(F,[tM ti],yM,op,aa); 
   warning('on','all');
   Yma  = y(end,:);  end

ymi = -Yma(CN);                   