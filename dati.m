% print DateTime

function dt = dati(tx)
c  = clock;  
dt = sprintf('%02d.%02d.%4d  %02d-%02d',c([3 2 1 4 5 ]));
fprintf('%50s\n%50s\n',dt,tx)