function [res] = f1(x)

    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    
    res = 2*(x1+x2+x3-3)^2 + (x1-x2)^2 + (x2-x3)^2;

end