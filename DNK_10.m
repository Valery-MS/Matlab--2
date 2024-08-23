% Search for periodic & almoust periodic(AP) oscillations in a DNA molecule

function DNK_10
global  a  b  y000  TE 
dati(1);
tic                     % Choice of methods for solving ODE
         % The methods № of: 1)dsm, m=6,8,10, 2)dop853, 3)ode113, 4)odex   
WIN = 2; % Whole Interval method Number(many points) of 1-3, best ode113
RPN = 3; % Refine Period method Number of 2-3, best ode113 (odex is deleted)
         % Odex has been removed: it does not work on (a, b), a>b (in fminbnd) 
         
aa = [0.0344 0.0446 0.0092 0.0144]; % a1 a2 d1 d2
CN = 1;                             % Coordinate Number of changed Initial Value 
y000 = zeros(1,4);                  % jig for initial values
a    = aa(1)/aa(3);   b = aa(2)/aa(4);  a2 = a*a;
if CN == 2,   o = a; a = b;   b = o;  end
opi  = 180/pi;                      % in radians
RelT = eps;  AbsT = eps; 
ht   = 1; %[1 0.5 0.25 0.1];        % t-step size

%                  Equilibrium positions Тетрадь15-77
%L  = 10*pi;  df = L/1000;    f  = 0:df:L;    sf = sin(f); 
%Q1 = a*sf-sin(f+asin(a*sf/b)); 
%2 = b*sf-sin(f+asin(b*sf/a));   plot(f,[Q1,Q2]);

if WIN == 1        % Whole Interval method Number=1 => WImet=dsm (default ds10)
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
   WIop  = DSet( ht, RelT,AbsT,Lmn, zam);
else
   if     WIN == 2,  WImet = @dop853; WIset = @dopset; 
   elseif WIN == 3,  WImet = @ode113; WIset = @odeset; 
   else,  errordlg('Wrong value of a WImet');end    
   WIop = WIset('RelTol',RelT,'AbsTol',AbsT,'InitialStep',ht);end 

if     RPN == 2,  RPmet = @dop853; RPset = @dopset; 
elseif RPN == 3,  RPmet = @ode113; RPset = @odeset; 
else,  errordlg('Wrong value of a RPMet');end  

RPop= RPset('RelTol',RelT,'AbsTol',AbsT,'InitialStep',[]);
erA = @(u,v) max(abs(u-v));
erR = @(u,v) max(abs(u./v-1));  
foE = @(f,a,b) a*sin(f)-sin( f+asin(a*sin(f)/b) );
ep2 = eps*2;  
t0  = 0;    t9 = 6000; 
t   = (t0:ht:t9)';          
MI  = 25;
RN  = 200;                     % root number
Df  = 1e-4;   Df2 = Df/2;      % Distance between roots f
Tab = nan(RN,7);
nh  = 4;                       Tabs = cell(nh,3);
hma = a+b-1;                   hmi  = a+b-0.5*(a/b+ a*b+b/a);  h0 = 1;
dh_ = (hma-h0)/nh;             htol = 0.01*dh_;
dh  = (hma-h0-2*htol)/nh;      h    = h0+htol-dh;
e1  = 1e-7;   e2 = 1e-5;       e3   = 1e-7; 
inf = sprintf('%s %s CN=%d t9=%d hma=%5.3g er=%6.1g %6.1g %6.1g',...
      func2str(WImet),func2str(RPmet),...
      CN,t9,hma,e1,e2,e3);  fprintf('%s\n',inf);  
wid = 3;                  % width of (A,B)-interval split ~ 3..5 in degrees
                          % Main loop
for i = 1:nh              % h-counter
  tic
  h  = h+dh;                    
  rD = sqrt(2*a*b*(h-hmi));     xc  = a*(h-b)+b;  
  xmi= (xc-rD)/a2;              xma = (xc+rD)/a2; 
  A0 = acos(1-max(0,xmi))+ep2;  B0  = acos(1-min(1,xma))-ep2;   BA0 = B0-A0;
  hQ = BA0/round(BA0*opi/wid);  Q   = B0:-hQ:A0+hQ;           % Q - queue
  A  = A0;                      B   = A0+hQ/2;
  fprintf('------ i=%d  h =%5.3g   %6.3g < fi10 <%5.3g\n',i,h,[A0 B0]*opi);
  c = 0;                      % root counter
  key=1;
  while A < B0-Df
     c = c+1; 
     [f,dy,eF,O] = fminbnd( @(f) APSr(key,t,f,h,CN,aa,WImet,RPmet,WIop,RPop),...
     A,B,optimset('TolX',eps,'MaxIter',MI,'Display','off'));  
     if eF < 0, waitfor(errordlg(O.message));end
     % ТЕ is sought from APSr 
     fifi(f,h,CN);                                        warning('off','all');
     [t_,y,n_] = WImet(@F_DNK,[0 TE/2 TE],y000,WIop,aa);  warning('on','all');
     dc = max(abs(y(2,3:4)));  dT = erA(y(3,1:2),y000(1:2));                                                    
     YE = dy<e1 && dc<e2 && dT<e3;
     Tab(c,:) = [ f*opi TE dy dc dT O.iterations YE]; 
     fprintf('%d %7.3f %8.1f %6.1g %6.1g %6.1g  %d %d\n',c,Tab(c,:));
     if YE && f-A > Df && B-f > Df
        B = f-Df;  Q = [Q f];
     else
        A = Q(end)+Df; Q(end)=[]; if isempty(Q),break;end
        B = Q(end)-Df; end,end
     
  Tab = sortrows( Tab(1:c,:), 1);
  Tabs(i,:) = {Tab inf sprintf('h=%5.2g  nr=%d  time=%dm\n',h,c,round(toc/60))}; 
  PRez( Tabs(i,:),7,e1);end 
toc
if exist('fig','var'), delete(fig); end
save(['DNK_10 ' date],'Tabs');