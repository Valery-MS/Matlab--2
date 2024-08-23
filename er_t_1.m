% Period error depending on t at CN coordinate of solution
% CN - Chaning Initial Value Index
% (tM,yM) - extremum point in discret solution y(:,CN)

function er = er_t_1(ti,iv,err,met,F,tM,yM,CN,op,varargin)

if ti == tM, er = err(yM(CN),iv);
else      
   op.InitialStep = (ti-tM)*0.5;
   warning('off','all');
   [t_, y] = met(F,[tM ti],yM,op,varargin{:});
   warning('on','all');
   er  = err(y(end,CN),iv);  end                % Abs or Relative err