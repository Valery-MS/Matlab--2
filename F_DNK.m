% Right hand function of DNK equation
% применяется в паре с функцией DNK_M в методах DS31n, DS11, DS21, DS31
% f = F(y)
% This variant is more than 2 times faster than with preliminary 
% allocation  f=zeros(4,1)

function F = F_DNK(t,y,a)

s = sin(y(1)+y(2));
F = [ y(3); y(4); -a(1)*sin(y(1))+a(3)*s; -a(2)*sin(y(2))+a(4)*s ];
