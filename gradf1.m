function [gf] = gradf1(x)

gf = 4*(x(1)+x(2)+x(3)-3)*ones(3,1) +[ 2*(x(1)-x(2)) ; -2*x(1)+4*x(2)-2*x(3) ; -2*(x(2)-x(3))];

end 