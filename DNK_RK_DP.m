% Solving of DNK equation by Runge-Kutta or Dormand-Prince methods
function DNK_RK_DP(ak,tt,metset,met,handles)
global a y0 t0 tn h kt RelT AbsT  
global tt y12 ak

op = metset('RelTol',RelT,'AbsTol',AbsT,'InitialStep',h);
warning('off','all')
tic; 
d = y0(1);
y01s = d:d:2*pi;
for y01 = y01s
   y0(1) = y01;
   [t, y] = met(@DNKf,tt,y0,op,a);
   y12 = y(:,1:2); 
   figure
   den = round(pi/y01);
   zer = zeros(size(t));
   subplot(3,1,1);  plot(t,y12(:,1),'b',t,zer,'k'); 
   meth = handles.Methods.String{handles.Methods.Value};
   dkm  = sprintf('%d, kt=%d  %.0s',den,kt,meth);
   title({['\phi_{10}=\pi/' dkm] '\phi_1'})
   subplot(3,1,2);  plot(t,y12(:,2),'r',t,zer,'k');    title('\phi_2')
   subplot(3,1,3);  plot(t,y12, t,zer,'k'); title('\phi_1, \phi_2')
   xlabel('t'); 
   title('\phi_1, \phi_2');  end
   
T = toc;
warning('on','all'); 
save(handles.DNK.UserData{1},'a','y0','t0','tn','h','kt','RelT','AbsT');