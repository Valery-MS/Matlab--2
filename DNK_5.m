% Search for periodic & almoust periodic(AP) oscillations in a DNA molecule

function DNK_5
global yTA                    % calcs in Tj = fminbnd(@(t) fnc(t,...)
global Figs                   % for fast deleting of figures  
global J bn TE ert Ts TG TT   % грубое точное
J  = 0;  Ts = nan(1000,1);  TT = nan(2,1000); TG = TT; ert = TT; 

tic                           %  Source data
SZ  = [459 230 926 515]; % = X Y Width Height 
Ffi = 1; %Flag fi 1)CN 2)1,2 3)1,2,3,4 4)1+2,3+4 5)1,2 vs Lin 6)1 vs Lin
qBs = [5 2 1 2 2 1 4 1];  % Qty of Blocks in Figure for Ffi (array)
qGs = [1 2 4 2 2 4 1 2];  % Qty of Grafics in Block for Ffi (array)
GRAF1 = 0;                % intermediate grafs are plottinivg or no (1/0)
GRAF2 = 0;                % 2nd part of grafs are plotting or no  (1/0)

a = [ 0.0344    0.0446    0.0092    0.0144];
t0  = 0;    tf = 6000;   
tA  = 3000; tB = 6000;   tik = 500;
T2ma = 100;     % if T2>T2ma, grafic is plotted. If T2ma->inf, T-> period
mp_t = 1;       % multiplier of t
        
                     % Choice of methods for solving ODE
       % The methods № of: 1)dsm, m=6,8,10, 2)dop853, 3)ode113, 4)odex   
WIN = 2; % Whole Interval method Number(many points) of 1-3, best ode113
RPN = 3; % Refine Period method Number of 2-3, best ode113 (odex is deleted)
% NOTE: Odex has been removed because it does not work on the interval (a, b),
% where a>b, which is used in @-func in fminbnd

y0 = zeros(1,4);   % is a dot of initial values
CN = 1;            % Coordinate Number where Initial Value is changed - IV(1)
if CN <= 2
   ivd0 = 1;         hivd = 0.1;       ivd9 = 50;    pi1 = pi/180;
   iv0  = ivd0*pi1;  hiv  = hivd*pi1;  iv9  = ivd9*pi1;
else
   iv0  = 0.01;        hiv  = 0.01;        iv9  = 1; end
niv = round((ivd9-ivd0)/hivd);
ivT = nan(niv+100,1);

              % initial deviation range - 1st pendulum: CN = 1
%ivd = ...;   % initial deviation range - 2st pendulum: CN = 2
%iv = ...;    % initial velocity  range - 1st pendulum: CN = 3
%iv = ...;    % initial velocity  range - 2st pendulum: CN = 4

aiv0 = abs(iv0);  if aiv0 == 0, aiv0 = 1; end

RelT = eps;  AbsT = aiv0*eps; 
Tol = 2.5e-3;                 % max Relative Tolerance of solution y on pattern
h   = 1; %[1 0.5 0.25 0.1];   % t-step size
eRA = true;                   % Rel or Abs err of y = ivi
if eRA, err = @(u,v) abs(u./v-1);   
else,   err = @(u,v) abs(u-v); end 

%*************************************************************************
Figs = [];
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
   else   errordlg('Wrong value of a WImet');end    
   WIop = WIset('RelTol',RelT,'AbsTol',AbsT,'InitialStep',h);end 

if     RPN == 2,  RPmet = @dop853; RPset = @dopset; 
elseif RPN == 3,  RPmet = @ode113; RPset = @odeset; 
else   errordlg('Wrong value of a RPMet');end  

RPop = RPset('RelTol',RelT,'AbsTol',AbsT,'InitialStep',[]);
Mets = sprintf('%s %s h=%g GRAF1,2=%d %d',...
        func2str(WImet),func2str(RPmet),h,GRAF1,GRAF2);

%*************************************************************************
t   = (t0:h:tf)';  % t0f = [t0 tf];
nt  = numel(t); 
%tP = [t; t(end)]; % used instead of t if t(end) is a multiple of the period 
t_  = t/mp_t;  XLi = [tA tB];  t9_= tf/mp_t;   XLi_= XLi/mp_t;  tik_= tik/mp_t;
zer = zeros(size(t));  % n = numel(t);     

                             % Source data output 
IB = round((tA-t0)/h)+1:round((tB-t0)/h)+1;  tIB = t(IB);

if niv == 1, qB = 1;
else         qB = qBs(Ffi); end  % qty of Grafics  for y plot

qG    = qGs(Ffi);       % qty of Grafics  for y plot
qGB   = qG*qB;
fi    = {'$\varphi_1$',    '$\varphi_2$',    '$v_1$',    '$v_2$'};
fiCNs = {'$\varphi_{10}$', '$\varphi_{20}$', '$v_{10}$', '$v_{20}$'};
fiCN  = fiCNs{CN};
                     
T   = Inf(2,niv);   % 1st 2 almoust periods(AP) & their difference for each iv
ivT = nan(1,niv);  ivTd = ivT; TEs = ivT;
ivL = iv0-hiv;  iv = iv0; 
i   = 0;  j = 1;
dt  = (t(2)-t(1))*0.7;
dn  = 10; 
MI  = 22;
wrn = false;       % on/off warnning message
erMa = 1e-10;

fprintf(...
'%s\nFfi=%d CN=%d ivd=%.2f:%.2f:%d Tol=%g tf=%d tA=%d tB=%d T2ma=%d dn=%d\n',...
Mets, Ffi, CN, ivd0,hivd,ivd9,Tol,tf,tA,tB,T2ma,dn);  

%                     Search for almost periods
while true  
    i = i+1;  ivL = ivL+hiv;  iv = iv+hiv;    
    [y, T(:,i)] = APSz(t,iv,err,y0,CN,a,Tol,WImet,WIop); 
    if all( T(:,i) == Inf )
       fprintf('i=%d ivd=%3.1f was skipped i.e n_=empty => T=[Inf Inf]\n',...
                i, iv/pi1);   
       continue;  end
    
    bn = 0;                       % 0;1 or t-boundaries numbers in APSr         
    [ivT(j) dyT exF ou] = fminbnd...
            (@(iV) er_iv(t,iV,err,y0,CN,dt,dn,a,WImet,RPmet,WIop,RPop),...
            ivL,iv,optimset('TolX',2*eps,'MaxIter',MI,'Display','off'));  %   
    if dyT > erMa || exF == 0
       if wrn, fprintf('skip %d  %d  %6.1f  %5.1g   %d\n',...
               i, j, iv/pi1, dyT, ou.iterations); end
       continue, end
    
    TEs(j) = TE;         % TE(1) - exact T for init. val. ivT(j)
    qTE = floor(t(end)/TE(1));
    tTE = [t(1) TE(1)*(1:qTE)];
    y0(CN) = ivT(j);
    warning('off','all');
    [tE,yT,nuf] = WImet(@RHFN,tTE,y0,WIop,a,4);        
    warning('on','all');
                    % Search of numbers n: tTE(n)=TE, 2TE, 3TE,...qTE*TE
    Ns = zeros(1,qTE);   
    k  = 1;
    for r=1:numel(tE), if abs(tE(r)-k*TE(1))<1e-10, Ns(k)=r; k=k+1;end, end

    ivTd(j) = ivT(j)/pi1;
    dyTs = max(abs(yT(Ns,CN)-ivT(j)));
    fprintf('%d  %d  %8.3f  %8.2f  %5.1g  %5.1g  %d\n',...
            i,j,ivTd(j),TE(1),dyT,dyTs,ou.iterations); 
    j = j+1; 
    if iv >= iv9, break,end,end  

ivd = ivT/pi1;
figure('Name',Mets,'NumberTit','on')
dim = size(T,2); 
plot(ivd(1:dim),T(1,:),'k', ivd(1:dim),T(2,:),'b', ivTd,TEs,'r.');
set(gca,'XGrid','on','YGrid','on','XMinorTick','on','TickDir','in');
ylabel('T_1, T_2'); xlabel(fiCN,'interpreter','latex')
fprintf('qiv=%d toc1=%g\n',niv,toc);

 TEs(:,1:j)

%save(['C:\Users\solop\Documents\MATLAB\ODE\COMMODE\' ...
%handles.DNK.UserData{1}],'SpecInpD','GenInpD');