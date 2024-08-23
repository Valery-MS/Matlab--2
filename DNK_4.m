% Search for periodic & almoust periodic(AP) oscillations in a DNA molecule

function DNK_4
global yTA   % calcs in solt с.164:  Tj = fminbnd(@(t) -metR(t,... )
global Figs  % for fast deleting of figures
global TE
tic 
                         %  Source data
%*************************************************************************
SZ  = [459 230 926 515]; % = X Y Width Height 

Ffi = 1; %Flag fi 1)CN 2)1,2 3)1,2,3,4 4)1+2,3+4 5)1,2 vs Lin 6)1 vs Lin
qBs = [5 2 1 2 2 1 4 1];  % Qty of Blocks in Figure for Ffi (array)
qGs = [1 2 4 2 2 4 1 2];  % Qty of Grafics in Block for Ffi (array)
GRAF1 = 0;                % intermediate grafs are plotting or no (1/0)
GRAF2 = 0;                % 2nd part of grafs are plotting or no  (1/0)

a = [ 0.0344    0.0446    0.0092    0.0144];

t0  = 0;    tf = 6000;   
tA  = 3000; tB = 6000;   tik = 500;
T2ma = 100;    % if T2>T2ma, grafic is plotted. If T2ma->inf, T-> period
mp_t = 1;       % multiplier of t

          
                     % Choice of methods for solving ODE
       % The methods № of: 1)dsm, m=6,8,10, 2)dop853, 3)ode113, 4)odex   
WIN = 2; % Whole Interval method Number(many points) of 1-3, best ode113
RPN = 3; % Refine Period method Number of 2-3, best ode113 (odex is deleted)
SPN = 3; % Second Period method Number of 1-3, best ode113  
% NOTE: Odex has been removed because it does not work on the interval (a, b),
% where a>b, which is used in solt in fminbnd

y0 = zeros(1,4);   % is a dot of initial values
CN = 1;      % Coordinate Number where Initial Value is changed - IV(1)
if CN <= 2
   ivd0 = 1;          hivd = 0.1;         ivd9 = 90;    pi180 = pi/180;
   iv0  = ivd0*pi180;  hiv  = hivd*pi180;  iv9  = ivd9*pi180;
else
   iv0  = 0.01;        hiv  = 0.01;        iv9  = 1; end
niv = round((ivd9-ivd0)/hivd)+1;
iv = nan(niv+100,1);

              % initial deviation range - 1st pendulum: CN = 1
%ivd = ...;   % initial deviation range - 2st pendulum: CN = 2
%iv = ...;    % initial velocity  range - 1st pendulum: CN = 3
%iv = ...;    % initial velocity  range - 2st pendulum: CN = 4

div = hiv/2;
aiv0 = abs(iv0);  if aiv0 == 0, aiv0 = 1; end

RelT = eps;  AbsT = aiv0*eps; 
Tol = 2.5e-3;                 % max Relative Tolerance of solution y on pattern
key = 3;                    % controls the kind of error in er_t, er_iv 
DT  = 2;                    % minimum period difference
kT  = 2;                    % the number of smallest periods
h   = 1; %[1 0.5 0.25 0.1]; % t-step size

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

%*************************************************************************
%if RPN == 4
%   RPmet = @odex; 
%   ITOL  = 0;     WORK  = zeros(1,9);   IWORK = zeros(1,71); 
%   IWORK(7) = 4;  IWORK(1) = 15000;    % NMAX, default = 1e4
%   RPop = struct('InitialStep',[],'RelTol',RelT,'AbsTol',AbsT,...
%                    'ITOL',ITOL, 'WORK',WORK, 'IWORK',IWORK);

if     RPN == 2,  RPmet = @dop853; RPset = @dopset; 
elseif RPN == 3,  RPmet = @ode113; RPset = @odeset; 
else   errordlg('Wrong value of a RPMet');end  

RPop = RPset('RelTol',RelT,'AbsTol',AbsT,'InitialStep',[]);

%*************************************************************************
if SPN == 1
   if SPN ~= WIN
      NEq  = 1; m = 3; Nh =1; ks = 1; kLmn = ks(1);
      switch NEq
        case 1, load('BestC_hs_Bess_InpD','BestC','hs');
        case 2, load('BestC_hs_Apry_InpD','BestC','hs');
        case 3, load('BestC_hs_Shrod_InpD','BestC','hs');
        otherwise, errordlg('Wrong value of a file number'); end  
      dsm   = { @ds6  @ds8  @ds10};    
      SPmet = dsm{m};
      zam   = [2 1 2]; % zam(1:2) iteration changes; zam(3)=1/2 if dop853/ode113
      Lmn   = BestC{m}{Nh}(kLmn,1:m);   
      SPop  = DSet( h, RelT,AbsT,Lmn, zam);
   else
      SPmet = WImet;
      SPop  = WIop;   end
else
   if     SPN == 2,  SPmet = @dop853; metset = @dopset; 
   elseif SPN == 3,  SPmet = @ode113; metset = @odeset; 
   else   errordlg('Wrong value of a WImet');end    
   SPop = metset('RelTol',RelT,'AbsTol',AbsT,'InitialStep',h);end 

Mets = sprintf('%s %s %s  h=%g GRAF1,2=%d %d',...
        func2str(WImet),func2str(RPmet),func2str(SPmet),h,GRAF1,GRAF2);


%*************************************************************************
t   = (t0:h:tf)';  % t0f = [t0 tf];
nt  = numel(t); 
%tP = [t; t(end)]; % used instead of t if t(end) is a multiple of the period 
t_  = t/mp_t;  XLi = [tA tB];  t9_= tf/mp_t;   XLi_= XLi/mp_t;  tik_= tik/mp_t;
MT  = 10:17;  % мн-ли периодов: различие TL и T, only for Lin Thai: Ffi>=5
zer = zeros(size(t));  % n = numel(t);     


                             % Source data output
fprintf(...
'%s\nFfi=%d CN=%d ivd=%.2f:%.2f:%d Tol=%g key=%d tf=%d tA=%d tB=%d T2ma=%d\n',...
Mets, Ffi, CN, ivd0,hivd,ivd9,Tol,key,tf,tA,tB,T2ma);

if Ffi >= 5             % for a Linear example of Thai
   d1=a(3); d2=a(4); e1=a(1)-d1; e2=a(2)-d2; ep=e1+e2; em=e1-e2;
   ed = sqrt(em^2+4*d1*d2);
   w1 = sqrt((ep+ed)*0.5);  w2 = sqrt((ep-ed)*0.5); 
   dw = w1-w2;  TL = 2*4*pi/dw; % tz = pi/dw;
   A = [[0 0 1 0]; [0 0 0 1]; [-e1 d1 0 0]; [d2 -e2 0 0]];
   [V, D] = eig(A);
   a_= real(V(1,1)); b_= real(V(2,1)); c_= imag(V(3,1)); d_= imag(V(4,1));
   e_= real(V(1,3)); f_= real(V(2,3)); g_= imag(V(3,3)); h_= imag(V(4,3));
   cow1 = cos(w1*t); cow2 = cos(w2*t);
   af = a_*f_; be = b_*e_; bf= b_*f_;   den = af-be;
   dwt = (w1-w2)*0.5*t;
   og_ = 2/den*[af*cos(dwt) bf*sin(dwt)];   % огибающая ф-ции yL_
   yL_ = [(af*cow1-be*cow2),  bf*(cow1-cow2)]/den;
   %yL = [cow1+cow2, -cow1+cow2]; - неправильно
   %cop = cos((w1+w2)*0.5*t);     - неправильно
   end
   
IB = round((tA-t0)/h)+1:round((tB-t0)/h)+1;  tIB = t(IB);

if niv == 1, qB = 1;
else         qB = qBs(Ffi); end  % qty of Grafics  for y plot

qG  = qGs(Ffi);       % qty of Grafics  for y plot
qGB = qG*qB;
wf  = round(78/h)+1;    % width of fractal 
nf  = round(20/h)+1:wf; % number of points on fractal
fi    = {'$\varphi_1$',    '$\varphi_2$',    '$v_1$',    '$v_2$'};
fiCNs = {'$\varphi_{10}$', '$\varphi_{20}$', '$v_{10}$', '$v_{20}$'};
fiCN  = fiCNs{CN};
                     % kT = 2
T   = Inf(kT,niv);   % 1st 2 almoust periods(AP) & their difference for each iv
i   = 1;             % ivd counter
ivi = iv0; 
YE  = false;         % YES: root perid intervcal is found
Bi  = nan(2,niv);    % Boundary initial values
hivi = hiv;
ivT  = nan(niv,1); 
ivTd = ivT;
TEs  = zeros(kT,niv);
j = 0;
cInf = 0;
izmh = 0;
               % Search for almost periods
while true   % figures loop
   if GRAF1  % intermediate grafs are plotting
      Figs = [Figs  figure('Positi',SZ,'Name',Mets,'NumberTit','on')];
      nF1 = Figs(end).Number; end
   if GRAF2
      Figs = [Figs  figure('Positi',SZ,'Name',Mets,'NumberTit','on')];
      nF2 = Figs(end).Number; end
  
   s = 1; 
   while s <= qB   % s-th plot Block in Fugure   
      [y, T(:,i)] = APS(ivi,t,key,y0,CN,h,nt,a,Tol,kT,...
                        WImet,RPmet,SPmet,WIop,RPop,SPop); 
                              
      if i > 1                            % SEARCH for bounderies of periods 
         TkTi = T(kT,i);
         if TkTi == Inf, cInf = cInf+1; end
         UP = TkTi - T(kT,i-1);
         if abs(T(1,i)-T(1,i-1)) < DT     % no jump to other period   
            if ~YE && UP < 0
               YE = true;   j = j+1;  
               if cInf
                  if i-cInf >= 1, ivL = iv(i-cInf-1); 
                  else            ivL = iv(1); end
                  cInf = 0;
               else     ivL = iv(i-2);    end
               Bi(:,j) = [ivL; ivi];
               
               [ivT(j) dyT]=fminbnd(@(iv) er_iv(t,iv,key,y0,CN,h,nt,a,Tol,kT,...
                   WImet,RPmet,SPmet,WIop,RPop,SPop),...
                   ivL,ivi,optimset('TolX',eps));
                   
               TEs(:,j) = TE;         % TE(1) - exact T for init. val. ivT(j)
               qTE = floor(t(end)/TE(1));
               tTE = [t(1) TE(1)*(1:qTE)];
               warning('off','all'); 
               y0(CN) = ivT(j);
               [tE,yT,nuf] = WImet(@RHFN,tTE,y0,WIop,a,4);        
               warning('on','all');
                    % Search of numbers n: tTE(n)=TE, 2TE, 3TE,...qTE*TE
               Ns = zeros(1,qTE);   
               k  = 1;
               for r = 1:numel(tE)
                  if abs(tE(r)-k*TE(1))<1e-10, Ns(k) = r; k = k+1;end,end

               ivTd(j) = ivT(j)/pi180;
               dyTs = max(abs(yT(Ns,CN)-ivT(j)));
               fprintf('j=%0.2d %7.3f %8.2f %5.1g %5.1g\n',j,ivTd(j),TE(1),dyT,dyTs);
         
               if GRAF2            % Graphs for init.value iv(i)                   
                  figure(nF2);   
                  ax = subplot(qGB,1,s); plot(tE,yT(:,CN),'k') %,t,zer,':k');
                  tTT = sprintf('%.5g $\\quad$',TE(~isnan(TE(1:kT))));  
                  title(['\textbf{\boldmath{' sprintf('%g  %s  %s', ...
                  ivTd(j),'$^\circ\quad$',tTT) '}}'],'interpreter','latex');
                  ylabel(fi{CN},'interpreter','latex');
                  set(ax,'XTick',0:tik_:t9_,'XLi',XLi_); 
                  s = s+1;  end
    
            elseif YE && UP > T(1,i)
               YE = false; end
         else  %if ~YE                     % jump, YE=0, search previous period   
            if hivi == hiv,   Ts = T(:,i); ivs = ivi; end      % s - save                 
            hivi = hivi/2;  ivi = ivi-hivi;  
            izmh = izmh+1;
            continue;       end,end
                                                    
      if GRAF1                  % intermediate plotting
         fi0 = sprintf('%g',iv(i)/pi180); 
         APPlot(s,i,nF1,qGB,t_,CN,y,Ffi,T,fi0); 
         s = s+1;  end
      
      iv(i) = ivi; i = i+1;  
      if   hivi == hiv,  ivi = ivi+hiv;
      else hivi =  hiv;  iv(i) = ivs; T(:,i) = Ts; i = i+1;  ivi = ivs+hiv; end

      if ivi > iv9, break,end,end  % while s<=qB - end of plotting of Block № s
  
   if GRAF1, xlabel('$t,[\times 1000],c$','interpreter','latex'); end
   if ivi > iv9, break,end,end  % whiletrue - end of plotting of Figure № nF1

ivd = iv/pi180;
figure('Name',Mets,'NumberTit','on')
dim = size(T,2); 
plot(ivd(1:dim),T(1,:),'k',ivd(1:dim),T(2,:),'b');
set(gca,'XGrid','on','YGrid','on','XMinorTick','on','TickDir','in');
ylabel('T_1, T_2'); xlabel(fiCN,'interpreter','latex')
fprintf('qiv=%d toc1=%g\n',niv,toc);

 TEs(:,1:j), Bi(:,1:j)

%save(['C:\Users\solop\Documents\MATLAB\ODE\COMMODE\' ...
%handles.DNK.UserData{1}],'SpecInpD','GenInpD');