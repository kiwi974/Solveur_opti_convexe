function [gf] = gradf2(x)

gf = [ -400*(x(1)*x(2)-x(1)^3) - 2*(1-x(1)) ; 200*(x(2)-x(1)^2)];

end 