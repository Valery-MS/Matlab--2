% Period error depending on t at CN coordinate of solution
% CN - Chaning Initial Value Index
% (tM,yM) - extremum point in discret solution y(:,CN)

function er = er_t(ti,met,F,tM,yM,CN,op,aa)
global Y

if ti == tM, er = -yM(CN);
else      
   op.InitialStep = (ti-tM)*0.5;
   warning('off','all');
   [t_,y] = met(F,[tM ti],yM,op,aa);
   warning('on','all');
   Y  =  y(end,:);
   er = -Y(CN);  end               