% 3 и 4 столбцы матрицы-функции DNK_M 

function f = DNK_34(t,y,varargin)

a = varargin{1};
s = sin(y(:,1)+y(:,2));
f = [-a(1)*sin(y(:,1))+a(3)*s, -a(2)*sin(y(:,2))+a(4)*s];