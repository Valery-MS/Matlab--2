% Search for periodic & almoust periodic(AP) oscillations in a DNA molecule

function DNK_6
global yTA                     % calcs in Tj = fminbnd(@(t) fnc(t,...)
global J bn TE ert Ts_ TG TT   % грубое точное

J  = 0;  Ts_ = nan(1000,1);  TT = nan(2,1000); TG = TT; ert = TT; 
                       
tic                     % Choice of methods for solving ODE
       % The methods № of: 1)dsm, m=6,8,10, 2)dop853, 3)ode113, 4)odex   
WIN = 2; % Whole Interval method Number(many points) of 1-3, best ode113
RPN = 3; % Refine Period method Number of 2-3, best ode113 (odex is deleted)
% NOTE: Odex has been removed because it does not work on the interval (a, b),
% where a>b, which is used in @-func in fminbnd

a  = [ 0.0344    0.0446    0.0092    0.0144];
y0 = zeros(1,4);   % is a dot of initial values
             % initial deviation range - 1st pendulum: CN = 1
%xd = ...;   % initial deviation range - 2st pendulum: CN = 2
%x = ...;    % initial velocity  range - 1st pendulum: CN = 3
%x = ...;    % initial velocity  range - 2st pendulum: CN = 4

CN = 1;            % Coordinate Number where Initial Value is changed - IV(1)
fiCNs = {'$\varphi_{10}$', '$\varphi_{20}$', '$v_{10}$', '$v_{20}$'};
fiCN  = fiCNs{CN};
if CN <= 2
   X0 = 1;       X9 = 90;       pi1 = pi/180;   % in degrees
   x0 = X0*pi1;  x9 = X9*pi1;   pi2 = 180/pi;   % in radians
else
   x0 = 0.01;    x9 = 1; end

ax0  = abs(x0);  if ax0 == 0, ax0 = 1; end
RelT = eps;  AbsT = ax0*eps; 
h    = 1; %[1 0.5 0.25 0.1];   % t-step size

if WIN == 1   % Whole Interval method Number=1 => WImet=dsm (default ds10)
   NEq  = 1; m = 3; Nh =1; ks = 1; kLmn = ks(1);
   switch NEq
      case 1, load('BestC_hs_Bess_InpD','BestC','hs');
      case 2, load('BestC_hs_Apry_InpD','BestC','hs');
      case 3, load('BestC_hs_Shrod_InpD','BestC','hs');
      otherwise, errordlg('Wrong value of a file number'); end      
   dsm   = { @ds6  @ds8  @ds10};    
   WImet = dsm{m};
   zam   = [2 1 2];  % zam(1:2) iteration changes; zam(3)=1/2 if dop853/ode113
   Lmn   = BestC{m}{Nh}(kLmn,1:m);   
   WIop  = DSet( h, RelT,AbsT,Lmn, zam);
else
   if     WIN == 2,  WImet = @dop853; WIset = @dopset; 
   elseif WIN == 3,  WImet = @ode113; WIset = @odeset; 
   else,  errordlg('Wrong value of a WImet');end    
   WIop = WIset('RelTol',RelT,'AbsTol',AbsT,'InitialStep',h);end 

if     RPN == 2,  RPmet = @dop853; RPset = @dopset; 
elseif RPN == 3,  RPmet = @ode113; RPset = @odeset; 
else,  errordlg('Wrong value of a RPMet');end  

RPop = RPset('RelTol',RelT,'AbsTol',AbsT,'InitialStep',[]);

ri   = 1e-3;                   % root indent - отступ от корня
Tol  = 2.5e-3;                 % max Relative Tolerance of solution y on pattern
eRA  = true;                   % Rel or Abs err of y = xi
if eRA, err = @(u,v) abs(u./v-1);   
else,   err = @(u,v) abs(u-v); end 

t0   = 0;    t9 = 6000; 
t    = (t0:h:t9)';  % t0f = [t0 t9];    
wrn1 = false;              % on/off warnning message
erMa = 1e-8; 
dt = (t(2)-t(1))*0.7;
NoJ = 5;                   % No Jump of period on NoJ steps
dn = 10; 
MI = 25;
RN = 200;                  % root number
A  = x0;  B = x9; 
Xs = nan(RN,1);            % roots x (init value) in degrees
Ts = Inf(RN,1);            % exact periods
j  = 0; 
Q  = B;                    % Q - queue
NP = 2;                    % Numb of Periods: periodicity check in T,2T...NP*T
Dx = 1e-3;                 % Distance between roots x
Cbn = 0;                   % Counter: how many times option bn=1 worked

inf = sprintf('%s %s CN=%d Tol=%g t9=%d NoJ=%d dn=%d',...
      func2str(WImet),func2str(RPmet),CN,Tol,t9,NoJ,dn);
fprintf('%s X=%d %d\n',inf,X0,X9);  

%                     Search for almost periods
while A < x9-ri                         
    for bn = 0:1                % 0;1 or t-boundaries numbers in APSr 
       bnf = bn;
       [x dy exF ou] = fminbnd...
            (@(X) er_iv(t,X,err,y0,CN,dt,NoJ,dn,a,WImet,RPmet,WIop,RPop),...
            A,B,optimset('TolX',eps,'MaxIter',MI,'Display','off'));  
       if dy < erMa
          if bnf, Cbn = Cbn+1; end
          break;  end,end
    
    xd = x*pi2;
    if dy > erMa || B-x < Dx || x-A < Dx
       Q(end) = [];  if isempty(Q), break; end 
       A = B;  B = Q(end);   
       if wrn1, fprintf('skip %d  %6.1f  %5.1g   %d\n',...
                j+1, xd, dy, ou.iterations); end
       continue, end
   
    j = j+1;
    y0(CN) = x;                                         warning('off','all');
    [tT,yT,n_]= WImet(@RHFN,[t0 (1:NP)*TE],y0,WIop,a,4); warning('on','all');
    dys = abs(yT(end-1:end,CN)-x);
    if dys(1) > erMa, S = sprintf('%d !!',exF); else, S = ''; end

    Xs(j) = x*pi2;  Ts(j) = TE;  B = x;  Q = [Q x];
    fprintf('%d  %8.3f  %8.2f  %5.1g  %5.1g  %5.1g  %d %s\n',...
    j,xd,TE,dy,dys,ou.iterations,S);  end  

XTs = sortrows( [ Xs(1:j) Ts(1:j) ], 1);
figure('Name',inf,'NumberTit','on')
plot( XTs(:,1), XTs(:,2), 'k.'); 
set(gca,'XGrid','on','YGrid','on','XMinorTick','on','TickDir','in');
ylabel('T'); xlabel(fiCN,'interpreter','latex')
title(sprintf('%d roots time=%gs Cbn=%d',j,toc,Cbn));
j