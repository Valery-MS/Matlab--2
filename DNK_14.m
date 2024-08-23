% Search for periodic & almoust periodic(AP) oscillations in a DNA molecule
% Comparison on period T ( by DY1 ) + fminsearch

function DNK_14
global  a  b  y000  TE 
dati
tic                     % Choice of methods for solving ODE
         % The methods № of: 1)dsm, m=6,8,10, 2)dop853, 3)ode113, 4)odex   
WIN = 2; % Whole Interval method Number(many points) of 1-3, best ode113
RPN = 3; % Refine Period method Number of 2-3, best ode113 (odex is deleted)
         % Odex has been removed: it does not work on (a, b), a>b (in fminbnd) 
         
aa  = [0.0344 0.0446 0.0092 0.0144]; % a1 a2 d1 d2
CN  = 1;                            % Coordinate Number of var Initial Value 
CNG = 3;                            % Coordinate Number of Gamma
y000 = zeros(1,4);                  % jig for initial values
a    = aa(1)/aa(3);   b = aa(2)/aa(4);  a2 = a*a;
if CN == 2,   o = a; a = b;   b = o;  end
opi  = 180/pi;                      % in radians
RelT = eps;  AbsT = eps; 
ht   = 1; %[1 0.5 0.25 0.1];        % t-step size

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
erA = @(u) max(abs(u));
erR = @(u,v) max(abs(u./v-1));   
 
t0  = 0;                      t9  = 6000; 
t   = (t0:ht:t9)';            MI  = 25;
RN  = 200;                    nh   = 10;     % root number    numb of h
Tab = nan(RN,7);              Tabs = cell(nh,3);             
hmi = a+b-0.5*(a/b+ a*b+b/a); hma  = a+b-1;          h0  = 1;
dh_ = (hma-h0)/nh;            htol = 0.01*dh_;       ep4 = eps*2; 
dh  = (hma-h0-2*htol)/nh;     h    = h0+htol-dh;
Df  = 1e-4;                   key  = 1;     % Distance between f & tolerance
e1  = 1e-4;   e2 = 1e-12;     EPS  = eps;             
inf = sprintf('%s %s CN=%d t9=%d hma=%5.3g er=%6.1g %6.1g EPS=%6.1g key=%d',...
      func2str(WImet),func2str(RPmet),...
      CN,t9,hma,e1,e2,EPS,key);  fprintf('%s\n',inf);  
                              % Main loop
for i = 1:nh                  % h-counter
  tic
  h  = h+dh;                    
  rD = sqrt(2*a*b*(h-hmi));     xc  = a*(h-b)+b;  
  xmi= (xc-rD)/a2;              xma = (xc+rD)/a2; 
  A  = acos(1-max(0,xmi))+ep4;  B   = acos(1-min(1,xma)); 
  Ao = A*opi;                   Bo  = B*opi; 
  fprintf('------ i=%d  h =%5.3g   %6.3g < fi10 <%5.3g\n',i,h,Ao,Bo);
  c = 0;  
  for f = [A+Df (ceil(Ao):floor(Bo))/opi B-Df]
     c  = c+1; 
     d0 = APSr(key,t,f,h,CN,aa,WImet,RPmet,WIop,RPop);   
     [fr,yr,eF,O] = fminsearch(@(F) DY1(F,h,CN,[t0 TE/2 TE],WImet,WIop,aa),f,...
     optimset('TolX',EPS));  
     fc = O.funcCount;   YE = d0<e1 && yr<e2;
     Tab(c,:) = [fr*opi TE d0 yr eF fc YE]; 
     fprintf('%d %7.3f %8.1f %6.1g %6.1g  %d  %d  %d\n',c,Tab(c,:));
     
     fifi(fr,h,CN);                                       warning('off','all');
     [t_,y,n_] = WImet(@F_DNK,[t0 TE/2 TE],y000,WIop,aa); warning( 'on','all');
     gc = y(2,3:4);          %erA(y(2,3:4));
     fT = y(3,1:2)-y000(1:2);  %erA(y(3,1:2)-y000(1:2));
     gT = y(3,3:4)-y000(3:4);  %erA(y(3,3:4)-y000(3:4)); 
     fprintf('%20.1g %6.1g   %6.1g %6.1g   %6.1g %6.1g\n',gc,fT,gT); end                                  

  Tab = sortrows( Tab(1:c,:), 1);
  Tabs(i,:) = {Tab inf sprintf('h=%7.4g  nr=%d  time=%dm\n',h,c,round(toc/60))}; 
  PRez( Tabs(i,:),4,1e-11);end 
toc
save(['DNK_14 ' dati],'Tabs');