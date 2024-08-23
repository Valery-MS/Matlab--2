% this is used in a case that should not be

function peril(tx1,J,tx2,t,y)
fprintf('%s=%6.2g %s :\n',tx1,J,tx2);
figure('Name','peril');
plot(t,y);