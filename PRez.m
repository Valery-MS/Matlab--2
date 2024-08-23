% Print of Rezults Tab, inf, 0, j, toc
function PRez(TI,col,er)

fiCN  = '$\varphi_{10}$';
Tab = TI{1};
figure('Name',TI{2}, 'NumberTit','on')
XT = Tab( Tab(:,col) < er, 1:2);
plot( XT(:,1), XT(:,2), 'k.'); 
set(gca,'XGrid','on','YGrid','on','XMinorTick','on','TickDir','in');
ylabel('T'); xlabel(fiCN,'interpreter','latex')
title(TI{3}) 