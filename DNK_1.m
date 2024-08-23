% Search for periodic & almoust periodic(AP) oscillations in a DNA molecule

function DNK_1
global yTA   % calcs in solt с.164:  Tj = fminbnd(@(t) -metR(t,... )
global Figs  % for fast deleting of figures

tic 
                         %  Source data
%*************************************************************************
SZ  = [459 230 926 515]; % = X Y Width Height 

Ffi = 5; %Flag fi 1)CN 2)1,2 3)1,2,3,4 4)1+2,3+4 5)1,2 vs Lin 6)1 vs Lin
qBs = [5 2 1 2 2 1 4 1];  % Qty of Blocks in Figure for Ffi (array)
qGs = [1 2 4 2 2 4 1 2];  % Qty of Grafics in Block for Ffi (array)
GRAF = 1;                 % transitional grafs are plotting or no

a = [ 0.0344    0.0446    0.0092    0.0144];

t0  = 0;    tf = 7000;   
tA  = 3000; tB = 7000;   tik = 500;
T2ma = 5000;    % if T2>T2ma, grafic is plotted. If T2ma->inf, T-> period
mp_t = 1;       % multiplier of t

          
                     % Choice of methods for solving ODE
          % The methods № of: 1)ds10, 2)dop853, 3)ode113, 4)odex   
WIN = 3;  % for Whole(biggest) Interval (1000s points): 1-3  best ode113
RPN = 4;  % to Refine the Period:                       2-4  best ode113,odex
SPN = 3;  % to search for the Second Period:            1-3  best ode113       

y0i = zeros(1,4);   % is a dot of initial values
CN  = 1;      % Coordinate Number where Initial Value is changed - IV(1)
ivd = [1:90]; % [5:45 48 50:90]; 
              % initial deviation range - 1st pendulum: CN = 1
%ivd = ...;   % initial deviation range - 2st pendulum: CN = 2
%iv = ...;    % initial velocity  range - 1st pendulum: CN = 3
%iv = ...;    % initial velocity  range - 2st pendulum: CN = 4
if CN <= 2, iv = ivd*pi/180; end       % range of fi-values(in radians)
iv1 = abs(iv(1));  if iv1 == 0, iv1 = 1; end

RelT = eps;  AbsT = iv1*eps; 
RT0 = 0.05;              % max Relative Tolerance of solution y on pattern
DT  = 2;                 % minimum period difference
kT  = 2;                 % the number of smallest periods
h   = 1; %[1 0.5 0.25 0.1]; 

%*************************************************************************
Figs = [];
if WIN == 1
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
if RPN == 4
   RPmet = @odex; 
   ITOL  = 0;     WORK  = zeros(1,9);   IWORK = zeros(1,71); 
   IWORK(7) = 4;  IWORK(1) = 15000;    % NMAX, default = 1e4
   RPop = struct('InitialStep',[],'RelTol',RelT,'AbsTol',AbsT,...
                    'ITOL',ITOL, 'WORK',WORK, 'IWORK',IWORK);
else
   if     RPN == 2,  RPmet = @dop853; RPset = @dopset; 
   elseif RPN == 3,  RPmet = @ode113; RPset = @odeset; 
   else   errordlg('Wrong value of a RfnMet');end    
   RPop = RPset('RelTol',RelT,'AbsTol',AbsT,'InitialStep',[]);end

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

Mets = sprintf('%s %s %s  h=%g',...
        func2str(WImet),func2str(RPmet),func2str(SPmet),h);


%*************************************************************************
qiv = numel(ivd);                      % qty of  initial values
t   = (t0:h:tf)';  % t0f = [t0 tf];
nt  = numel(t);
dt  = t(2)-t(1);   dt_5 = dt*0.5;  
%tP  = [t; t(end)]; % used instead of t if t(end) is a multiple of the period 
t_  = t/mp_t;  XLi = [tA tB];  t9_= tf/mp_t;   XLi_= XLi/mp_t;  tik_= tik/mp_t;
MT  = 10:17;   % мультипликаторы периодов: различие TL и T
zer = zeros(size(t));  % n = numel(t);     


                             % Source data output
fprintf('%s\nFfi=%d CN=%d ivd=%d:%d:%d RT0=%g tf=%d tA=%d tB=%d T2ma=%d\n',...
         Mets, Ffi, CN, ivd(1),ivd(2)-ivd(1),ivd(end),RT0,tf,tA,tB,T2ma);

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

if qiv == 1, qB = 1;
else         qB = qBs(Ffi); end  % qty of Grafics  for y plot

qG  = qGs(Ffi);       % qty of Grafics  for y plot
qGB = qG*qB;
wf  = round(78/h)+1;    % width of fractal 
nf  = round(20/h)+1:wf; % number of points on fractal
fi    = {'$\varphi_1$',    '$\varphi_2$',    '$v_1$',    '$v_2$'};
fiCNs = {'$\varphi_{10}$', '$\varphi_{20}$', '$v_{10}$', '$v_{20}$'};
fiCN  = fiCNs{CN};
                      % kT = 2
T  = Inf(2*kT-1,qiv); % 1st 2 almoust periods(AP) & their difference for each iv
Tm = nan(1,qiv);      % T for Min RTj
RT = T;               % Table of Relative Tolerances
i  = 0;               % ivd counter
                    
               % Search for almost periods
while true   % figures loop
   if GRAF
      Figs = [Figs  figure('Positi',SZ,'Name',Mets,'NumberTit','on')];
      nF = Figs(end).Number; end
   s = 1; 
   while s <= qB   % s-th plot Block in Fugure
%function        
      i   = i+1; 
      fi0 = sprintf('%g',ivd(i));
      ivi = iv(i);  y0i(CN) = ivi; 
      warning('off','all'); 
      [to,y,nuf] = WImet(@RHFN,t,y0i,WIop,a,4);        
      warning('on','all');          
      yi = y(:,CN);
      Nopp = find(ivi*yi<0,1); % Number of 1st element opposed to y0
      yi(1:Nopp) = NaN;        % Turn off value y0 from yi
      
                        % Double selection
      rt = abs(yi(Nopp+1:end)/ivi-1);  % rel tols of y values equal to the ivi
      n_ = find(rt<RT0);               % rt numbers: rt<RT0 (minimums ascending)
      rtn = rt(n_);                    % selection 1 by rt<RT0

      q = 0; x = 1; cs = 0;  I = n_;   % cs - cч соседей      
      nN = numel(n_);            % Исключаются соседние с минимумами точки -
      for j = 1:nN               % тот же min, но погрешность больше
         if j < nN  &&  n_(j+1)-n_(j) == 1       %  selection 2 
            if rtn(j+1) < rtn(j), x = j+1;  end  % by condition
            cs = cs+1;                           % neighbors(n_(j+1)-n_(j)=1)
            continue,end              
         if x < j-cs, x = j-cs; end              % are deleting
         q = q+1;  cs = 0;
         I(q) = n_(x);  end
      Ms = Nopp+I(1:q);     
           
      j = 0; r = 0; 
      while r < kT  &&  j < q    
         j = j+1;  M = Ms(j);   tM = t(M);
         if r > 0                  % search for a second (or more) AP
            dr = tM./T(1:r,i);
            if any(abs(dr-round(dr))<0.05), continue,end,end
        
         if     M == nt,           t1 = tM-dt_5;  t2 = tM+dt; 
         elseif yi(M-1) > yi(M+1), t1 = tM-dt_5;  t2 = tM;
         else                      t1 = tM;       t2 = tM+dt_5;    end
      
         % sign before dop853t = -sign(y0iC); more accurate AP value (Tj)
         % {'dop853' 'ode113_c' 'ode45_c' 'odex'}
         
         Tj = fminbnd(@(t) -solt(RPmet, @F_DNK,tM,t,y(M,:),CN,RPop,a),...
               t1,t2,optimset('TolX',eps));  
         warning('off','all');                          % for dop853, ode113
         Tfj = Tj :h: Tj+(wf-1)*h;                      % T-fractal
         [Tfj,yf] = SPmet( @F_DNK, Tfj, yTA, SPop, a); 
         warning('on','all');                           % for dop853, ode113
         RTj = max(abs(y(nf,CN)-yf(nf,CN)))/ivi;
         if RTj<RT0, r=r+1;  T(r,i)=Tj;  RT(r,i)=RTj; end,end
     
      Tm(i) = T(1,i);  
      if RT(2,i)<RT(1,i), Tm(i) = T(2,i); end
      T(kT+1:end,i) = T(2:kT,i)-T(1:kT-1,i);   
      
                              % Graphs for init.value iv(i)
      if GRAF &&  kT >= 2  &&  T(kT,i) > T2ma 
        figure(nF);    
        if Ffi == 1         % qGB=5*1, qG=1,  plot  fi(CN)        
          ax = subplot(qGB,1,s); plot(t_,y(:,CN),'k') %,t,zer,':k');
          tTT = sprintf('%.5g $\\quad$',T(~isnan(T(1:kT,i)),i));  
          %title(['\textbf{\boldmath{' sprintf('%s  %s',...
          %[fi0 '$^\circ\quad$'],tTT) '}}'],'interpreter','latex'); 
          %title(['\textbf{\boldmath{' sprintf('%s T = %s',[fiCN '=' fi0 ...
          title(['\textbf{\boldmath{' sprintf('%s  %s',[fi0 ...
          '$^\circ\quad$'],tTT) '}}'],'interpreter','latex');
          ylabel(fi{CN},'interpreter','latex');
          set(ax,'XTick',0:tik_:t9_,'XLi',XLi_);
        elseif Ffi == 2  
          for k = 1:qG    % qGB=2*2, plot fi1, fi2 
            ax = subplot(qGB,1,2*(s-1)+k); plot(t,y(:,k),'b',t,zer,'k'); 
            ylabel(fi{k},'interpreter','latex'); 
            set(ax,'XMinorGrid','on','XMinorTick','on','XTick',0:tik:tf);
            if k==1
               tTT = sprintf('%.5g $\\quad$',T(~isnan(T(1:kT,i)),i));
               title(['\textbf{\boldmath{' sprintf('%s T = %s',...
               [fiCN '=' fi0 '$^\circ\quad$'],tTT) '}}'],...
               'interpreter','latex');  end,end
        elseif Ffi == 3
          for k = 1:qG    % qGB=4*1, plot fi1,...fi4
            ax = subplot(qGB,1,k); 
            plot(t,y(:,k),'b',t,zer,'k'); 
            ylabel(fi{k},'interpreter','latex'); 
            set(ax,'XMinorGrid','on','XMinorTick','on','XTick',0:tik:tf);
            if k == 1
               tTT = sprintf('%.5g  ',T(~isnan(T(1:kT,i)),i)); 
               title(sprintf('%s (%s)',[fiCN '=' fi0 '$\quad$'],tTT),...
               'interpreter','latex');  end,end
        elseif Ffi == 4   % qG=2, plot fi1 U fi2, f3 U fi4
          for k = 1:qG  
            ax = subplot(qGB,1,2*(s-1)+k); k2 = k+k;
            plot(tIB,y(IB,k2-1),'b',tIB,y(IB,k2),'k',tIB,zer(IB),':k'); 
            ylabel([fi{k2-1} ',' fi{k2}],'interpreter','latex');
            set(ax,'XMinorGrid','on','XMinorTick','on','XTick',0:tik:tf,...
            'XLi',XLi_ );
            if k == 1
               tTT = sprintf('%.5g $\\quad$',T(~isnan(T(1:kT,i)),i));    
               title(['\textbf{\boldmath{' sprintf('%s  T = %s',...
               [fiCN '=' fi0 '$^\circ\quad$'],tTT) '}}'],...
               'interpreter','latex'); end,end       
        elseif Ffi == 5 % f1+L1, f2+L2 Comparison with linear approx (Thai)
          for k = 1:qG 
            yL = ivi*yL_; 
            ax = subplot(qGB,1,2*(s-1)+k); 
            plot(tIB,[y(IB,k),yL(IB,k)],tIB,zer(IB),':k',...
                [MT*TL MT*Tm(i)],0,'.k'); 
            ylabel([fi{k} ',' fi{k} 'Lin'],'interpreter','latex'); 
            set(ax,'XMinorGrid','on','XMinorTick','on','XTick',0:tik:tf,...
            'XLi',XLi_ ); 
            if k == 1
               tTT = sprintf('%.5g $\\quad$',T(~isnan(T(1:kT,i)),i)); 
               title(['\textbf{\boldmath{' sprintf('%s  T = %s',...
               [fiCN '=' fi0 '$^\circ\quad$'],tTT) '}}'],...
               'interpreter','latex');  end,end
        elseif Ffi == 6 % f1 vs L2 Comparison with linear approx (Thai)
          yL = ivi*yL_;
          tz=0;
          ax(1) = subplot(qGB,1,2*(s-1)+1); 
          T1 = T(1,i);  qT = round(tB/T1);
          plot(tIB,y(IB,CN),'b',tIB,zer(IB),':k',...
            tz+T1*(0:qT),0,'.k',tz+TL*(0:qT),0,'.r');
          ylabel(fi{CN},'interpreter','latex'); 
          tTT = sprintf('%.5g $\\quad$',T(~isnan(T(1:kT,i)),i));        
          title(['\textbf{\boldmath{' sprintf('%s  T = %s',...
          [fiCN '=' fi0 '$^\circ\quad$'],tTT) '}}'],'interpreter','latex'); 
        
          ax(2) = subplot(qGB,1,2*(s-1)+2); 
          plot(tIB,yL(IB,CN),'b',tIB,zer(IB),':k',...
             tz+T1*(0:qT),0,'.k',tz+TL*(0:qT),0,'.r'); 
          title(sprintf('\\omega_1=%g \\omega_2=%g \\TL=%g',w1,w2,TL))
          ylabel('$A(\cos w_1+\cos w_2$)','interpreter','latex');   
        elseif Ffi == 7 % |f1-L1| 
          yL = ivi*yL_; 
          ax = subplot(qGB,1,s); 
          plot(tIB,(y(IB,CN)-yL(IB,CN))/ivi) %,t,zer,':k');
          tTT = sprintf('%.5g $\\quad$',T(~isnan(T(1:kT,i)),i));  
          %title(['\textbf{\boldmath{' sprintf('%s  %s',...
          %[fi0 '$^\circ\quad$'],tTT) '}}'],'interpreter','latex'); 
          %title(['\textbf{\boldmath{' sprintf('%s T = %s',[fiCN '=' fi0 ...
          title(['\textbf{\boldmath{' sprintf('%s  %s',[fiCN '=' fi0 ...
          '$^\circ\quad$'],tTT) '}}'],'interpreter','latex');
          ylabel(fi{CN},'interpreter','latex');
          set(ax,'XTick',0:tik_:t9_,'XLi',XLi_);
        elseif Ffi == 8 % fi1, fi2, L1, L2, og        
          for k = 1:qG 
             yL = ivi*yL_; og = ivi*og_;
             ax = subplot(qGB,1,2*(s-1)+k); 
             plot(tIB,[y(IB,k),yL(IB,k),og(IB,k)],tIB,zer(IB),':k',...
                 [MT*TL MT*Tm(i)],0,'.k'); 
             ylabel([fi{k} ',' fi{k} 'Lin'],'interpreter','latex'); 
             set(ax,'XMinorGrid','on','XMinorTick','on','XTick',0:tik:tf,...
             'XLi',XLi_ ); 
             if k == 1
               tTT = sprintf('%.5g $\\quad$',T(~isnan(T(1:kT,i)),i)); 
               title(['\textbf{\boldmath{' sprintf('%s  T = %s',...
               [fiCN '=' fi0 '$^\circ\quad$'],tTT) '}}'],...
               'interpreter','latex');  end,end,end

        s = s+1; end % if GRAF=1 - end of plotting for i.v. iv(i) 
    
      if i == qiv, break,end,end  % while s<=qB - end of plotting of Block № s
  
   %xlabel('\textbf{\boldmath{$t,\times 10^3$}}','interpreter','latex');
   xlabel('$t,[\times 1000],c$','interpreter','latex') %xlabel('t,');
   
   if i == qiv, break,end,end  % whiletrue - end of plotting of Figure № nF


%{
figure('Name',sprintf('RelTol=%g',RT0))    %,'NumberTitle','off')
%ax1=subplot(2,1,1);
plot(ivd,T(1,:));
set(gca,'XGrid','on','YGrid','on','XMinorTick','on','TickDir','out');
%title(sprintf('1-й период  RT0=%g',RT0));
ylabel('T1'); xlabel(fiCN,'interpreter','latex');
%ax2=subplot(2,1,2);
figure
plot(ivd,T(2,:));
set(gca,'XGrid','on','YGrid','on','XMinorTick','on','TickDir','out');
%title('2-й период');
ylabel('T2');xlabel(fiCN,'interpreter','latex');
%ax3=subplot(3,1,3); plot(ivd,Tm);     title('комб период');ylabel('Tm')
%}

TAll = toc; qiv_T1_TAll = [qiv, TAll/qiv, TAll]

figure('Name',Mets,'NumberTit','on')
plot(ivd,T(1,:),'k',ivd,T(2,:),'b');
set(gca,'XGrid','on','YGrid','on','XMinorTick','on','TickDir','in');
ylabel('T_1, T_2'); xlabel(fiCN,'interpreter','latex')


                     % SEARCH for periods
nT  = find( T(kT,:)>T2ma );
qnT = numel(nT);
TLV = T(1:kT,nT);    % Table of all (Little & Very Different)AP, which T2>T2ma
nV  = nan(1,qnT); 
qV  = 0;             % qty of VD AP

                     % 1. Search for Very Different Almoist Periods (VD AP)
Bi = nan(2,qnT);     % Boundary initial values
cs = 0;              % neighbour counter
for j = 1:qnT        % Deleting of Little Different AP: TLV-TL=TV
   if j < qnT  &&  abs(TLV(1,j+1)-TLV(1,j)) < DT
      if TLV(kT,j+1) > TLV(kT,j), x = j+1; end
      cs = cs+1;
      continue
   else
      if x < j-cs, x = j-cs; end     
      qV = qV+1;   nV(qV) = x;
      if cs
         if     x == j,    Bi(1,qV) = iv(x-1); 
         elseif x == j-cs; Bi(2,qV) = iv(x+1); 
         else              Bi(1,qV) = iv(x-1); Bi(2,qV) = iv(x+1); end,end
      cs = 0;    
   end,end

                     % 2. Search for Boundaries of VD AP
nV   = nV(1:qV);
TB   = nan(2*kT-1,qV);
ivdV = ivd(nV);  
divs = [-div div]; 
ivT  = nan(1:qV);

for j = 1:qV 
   i   = nV(j);
   ivi = iv(i);
   T1  = T(1,i);
   for k = 1:2
      di = divs(k);  
      if isnan(Bi(k,j))
         while true
            ivB = ivi+di;  
            [y, TB(:,j), TBm] = APS(ivB,t,y0,CN,h,nt,a,RT0,kT,wf,nf,...
                                    WImet,RPmet,SPmet,WIop,RPop,SPop);  
            if TB(1,j)-T1 > DT,  di = 0.5*di; 
            else            Bi(1,j) = ivB; break; end,end,end,end  
 
   ivT(j) = fminbnd(@(iv) solD(iv,t,y0,CN,h,nt,a,RT0,kT,wf,nf,...
                    WImet,RPmet,SPmet,WIop,RPop,SPop),...
                    Bi(1,j),Bi(2,j),optimset('TolX',eps));  
end
%{

ns = ones(14,1);  r=1;  %  bounds of periods change  
for q=1:qiv-1
   if abs(T(1,q+1)-T(1,q))>10, r=r+1; ns(r)=q+1; end,end
ns(r+1)=qiv;
i=0; 
for nF=1:4
   figure
   s=0;
   while s<3  &&  i<r
      i=i+1; ii=ns(i):ns(i+1)-1;
      if ns(i)<ns(i+1)-1
         s=s+1; ax(s)=subplot(3,1,s); 
         plot(ii,T(1,ii)); 
         ylabel('T1');  end,end
   set(ax,'XGrid','on','YGrid','on','XMinorTick','on');end

nt1=numel(U1);     nt2=nt1+numel(U2); nt3=nt2+numel(U3); nt4=nt3+numel(U4);
nt5=nt4+numel(U5); nt6=nt5+numel(U6); nt7=nt6+numel(U7);  
N1=1:nt1; N2=nt1+1:nt2; N3=nt2+1:nt3; N4=nt3+1:nt4; N5=nt4+1:nt5;
N6=nt5+1:nt6; N7=nt6+1:nt7;
plot(U1,W(N1),U2,W(N2),U3,W(N3),U4,W(N4),U5,W(N5),U6,W(N6),U7,W(N7));
set(gca,'XGrid','on','YGrid','on','XMinorTick','on','TickDir','out');
ylabel('T1, T2'); xlabel(fiCN,'interpreter','latex')

%save(['C:\Users\solop\Documents\MATLAB\ODE\COMMODE\' ...
%handles.DNK.UserData{1}],'SpecInpD','GenInpD');
%}