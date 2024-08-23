% Получение nt0*dim-матрицы, используя DNK в nt0 точках t
% f_M = [ f(y(t1))'; f(y(t2))';...; f(y(t_nt0))' ], here ' is transposition
% применяется в паре с функцией DNK в методах DS31n, DS11, DS21, DS31
% y = nt0*dim - matrix,
% where nt0 = number of initial t-points for begin of solution calculating
%            ( at m=1,2,3 nt0 = 6, 12, 16 )             
%       dim = dimesion of system y' = F(y) 
%            ( here we consider dim = 2,4 )
% Transformation algorithm from F to Fv:
%  1) y(i) ->  y(:,i)
%  2) ";" -> ","  in formula of f
% Note: for nonautonom ODE t is n0-dimensional column,
%       => operatios '*' -> '.*', etc
% Here ODE is autonome

function F = F_DNK_M( t, y, a )

s = sin(y(:,1)+y(:,2));
F = [y(:,3), y(:,4), -a(1)*sin(y(:,1))+a(3)*s, -a(2)*sin(y(:,2))+a(4)*s];