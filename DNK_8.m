% Search for periodic & almoust periodic(AP) oscillations in a DNA molecule

function DNK_8
global  a  b  y000  TE 
dati(1)
tic                     % Choice of methods for solving ODE
         % The methods № of: 1)dsm, m=6,8,10, 2)dop853, 3)ode113, 4)odex   
WIN = 2; % Whole Interval method Number(many points) of 1-3, best ode113
RPN = 3; % Refine Period method Number of 2-3, best ode113 (odex is deleted)
         % Odex has been removed: it does not work on (a, b), a>b (in fminbnd) 
         
aa = [0.0344 0.0446 0.0092 0.0144]; % a1 a2 d1 d2
CN = 1;                             % Coordinate Number of changed Initial Value 
y000 = zeros(1,4);                  % jig for initial values
a    = aa(1)/aa(3);   b = aa(2)/aa(4);  a2 = a*a;
hmi  = a+b-0.5*(a/b+ a*b+b/a);  hma = a+b-1;
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

RPop = RPset('RelTol',RelT,'AbsTol',AbsT,'InitialStep',[]);
ri   = 1e-3;                   % root indent - отступ от корня
erA  = @(u,v) max(abs(u-v));
erR  = @(u,v) max(abs(u./v-1));   
  
t0 = 0;    t9 = 6000; 
t  = (t0:ht:t9)';  % t0f = [t0 t9];    
MI = 25;
RN  = 200;                    nh   = 4;     % root number    numb of h
Tab = nan(RN,7);              Tabs = cell(nh,3);             
hmi = a+b-0.5*(a/b+ a*b+b/a); hma  = a+b-1;          h0  = 1;
dh_ = (hma-h0)/nh;            htol = 0.01*dh_;       
dh  = (hma-h0-2*htol)/nh;     h    = h0+htol-dh;
Df  = 1e-3;                    % Distance between roots f
er1 = 1e-8;                   EPS  = eps;    ep2 = eps*2;
key = 1;        key2 = 0;   % choise of Root Interval (Rint) see p.57 Red Spine
inf = sprintf('%s %s CN=%d t9=%d hma=%5.3g er1=%6.1g EPS=%6.1g key=%d nh=%d',...
      func2str(WImet),func2str(RPmet),...
      CN,t9,hma,er1,EPS,key,nh);  fprintf('%s\n',inf);  

wid = 3;                      % width of (A,B)-interval split ~ 3..5 in degrees
k=0;                                % Main loop
for i = 1:nh                  % h-counter
  tic
  h  = h+dh;                   
  rD = sqrt(2*a*b*(h-hmi));     xc  = a*(h-b)+b;  
  xmi= (xc-rD)/a2;              xma = (xc+rD)/a2; 
  A  = acos(1-max(0,xmi))+ep2;  B0  = acos(1-min(1,xma))-ep2;   BA0 = B0-A;
  hQ = BA0/round(BA0*opi/wid);  Q   = B0:-hQ:A+hQ;            % Q - queue
  B  = Q(end);                  j   = 0;                      % root counter
  fprintf('------ i=%d  h =%5.3g   %6.3g < fi10 <%5.2g\n',i,h,[A B0]*opi);
  c = 0;                      % root counter
  while A < B0-ri 
    c = c+1; 
    [f,dy,eF,O] = fminbnd( @(f) APSr(key,t,f,h,CN,aa,WImet,RPmet,WIop,RPop),...
    A,B,optimset('TolX',EPS,'MaxIter',MI,'Display','off'));  
    fifi(f,h,CN);                                         warning('off','all');
    [t_,y,n_] = WImet(@RHFN,t0+[0 TE/2 TE],y000,WIop,aa,4);
                                                           warning('on','all');
    NOT = dy > er1 || B-f < Df || f-A < Df;
    if ~NOT && key2,k=k+1;end
    if NOT  || ~NOT && key2, key2 = mod(key2+1,3);end
    Tab(c,:) = [ ...
    f*opi TE dy max(abs(y(2,3:4))) erA(y(3,1:2),y000(1:2)) O.iterations ~NOT];
    fprintf('%d %7.3f %8.1f %6.1g %6.1g %6.1g   %d  %d\n',c,Tab(c,:));
    if     key2 == 1,          B = f;                   
    elseif key2 == 2,  A = B;  B = Q(end); Q(end) = [];    
    elseif NOT,                B = f;                   
    else,                      B = f;      Q   = [Q f]; end,end 

  Tab = sortrows( Tab(1:c,:), 1);
  PRez( Tab,h,inf,er1,c,toc);
  Tabs{i} = Tab;  end 
toc
if exist('fig','var'), delete(fig); end
save(['DNK_7 ' date],'Tabs','inf');
k
j