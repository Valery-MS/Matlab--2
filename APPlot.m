% Almoust Period Plotting

function APPlot(s,i,nF,qGB,t_,CN,y,Ffi,T,fi0)
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
      set(ax,'XMinorGrid','on','XMinorTick','on','XTick',0:tik:tf,'XLi',XLi_ );
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
      [MT*TL MT*T(nmi,i)],0,'.k'); 
      ylabel([fi{k} ',' fi{k} 'Lin'],'interpreter','latex'); 
      set(ax,'XMinorGrid','on','XMinorTick','on','XTick',0:tik:tf,'XLi',XLi_ ); 
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
      [MT*TL MT*T(nmi,i)],0,'.k'); 
      ylabel([fi{k} ',' fi{k} 'Lin'],'interpreter','latex'); 
      set(ax,'XMinorGrid','on','XMinorTick','on','XTick',0:tik:tf,'XLi',XLi_ ); 
      if k == 1
         tTT = sprintf('%.5g $\\quad$',T(~isnan(T(1:kT,i)),i)); 
         title(['\textbf{\boldmath{' sprintf('%s  T = %s',...
         [fiCN '=' fi0 '$^\circ\quad$'],tTT) '}}'],...
         'interpreter','latex');  end,end,end