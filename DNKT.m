% Search for periodic & almoust periodic(AP) oscillations in a DNA molecule
% for a given T, find y0

function DNKT
global  a  b  h y000  Y     
                       
tic                     % Choice of methods for solving ODE
         % The methods № of: 1)dsm, m=6,8,10, 2)dop853, 3)ode113, 4)odex   
WIN = 2; % Whole Interval method Number(many points) of 1-3, best ode113
RPN = 3; % Refine Period method Number of 2-3, best ode113 (odex is deleted)
         % Odex has been removed: it does not work on (a, b), a>b (in fminbnd) 
         
aa = [0.0344 0.0446 0.0092 0.0144]; % a1 a2 d1 d2
y000 = zeros(1,4);                  % jig for initial values
a    = aa(1)/aa(3);   b = aa(2)/aa(4);  a2 = a*a;
hmi  = a+b-0.5*(a/b+ a*b+b/a);  hma = a+b-1;
nh   = 5;
htol = 0.02*(hma-hmi)/nh;
dh   = (hma-hmi-2*htol)/nh;
%dh  = 0.1;                         % stepsize of energy integral
%nh  = floor((hma-hmi)/dh);
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
eRA  = 0;                           % Rel or Abs err of y = xi
if eRA, err = @(u,v) abs(1-u./v);   
else,   err = @(u,v) abs(  u -v); end 

t0  = 0;    t9 = 9000; 
t   = (t0:ht:t9)';            % t0f = [t0 t9];       
MI  = 25;
RN  = 10;                    % root number
Tab = nan(RN,5);
Tabs= cell(1,nh);
h   = 1+htol-dh;
inf = sprintf('%s %s t9=%d hma=%5.3g eRA=%d',...
      func2str(WImet),func2str(RPmet),t9,hma,eRA);  fprintf('%s\n',inf); 
fs = nan(1,RN);
hT = 1;        Ts = 0.5:hT:RN*hT; 
                              % Main loop
for i = 1:nh                  % h-counter
  tic
  h  = h+dh;                  
  rD = sqrt(2*a*b*(h-hmi));   xc  = a*(h-b)+b;  
  xmi= (xc-rD)/a2;            xma = (xc+rD)/a2; 
  A0 = acos(1-max(0,xmi))+eps;
  B0 = acos(1-min(1,xma))-2*eps;  
  j  = 0;                     % root counter
  fprintf('------ i=%d  h =%5.3g  %6.3g < fi10 <%6.3g\n',i,h,[A0 B0]*opi);
  
  for  T = Ts
    j = j+1; fprintf('T=%d\n',T);
    %AB = rint(fs,[t0 T],j,A0,B0,100,WImet,WIop,aa);
    %if isempty(AB),  continue, end
    [f,dy,exF,ou] = fminbnd( @(f) APST(t0,f,T,aa,WImet,WIop),...
    A0,B0,optimset('TolX',eps));      
    dys = max(err(Y(3,:),y000)); 
    if dys(end) < 1e-6
       if dys > 1e-11, eX = sprintf('%d !!',exF); else, eX = ''; end
       fd = f*opi;         
       Tabj = [ fd  TE max(abs(y(2,:))) dy  dys ];
       Tab(j,:) = Tabj; 
       fprintf('%d %7.3f %8.2f %7.1g %7.1g  %5.1g  %d %s\n',...
       j,Tabj,dys,ou.iterations,eX); end,end 

  Tab = sortrows( Tab(1:j,:), 1);
  PRez( Tab,inf,0,j, toc);
  Tabs{i} = Tab;  end 
toc
if exist('fig','var'), delete(fig); end
save(['DNK_7 ' date],'Tabs','inf');
j