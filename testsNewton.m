% Fichier de tests de l'algorithme de Newton, sans la méthode des 
% régions de confiance

res = {};

% **** Phase 1 ****
disp(['*************************************************']);
disp(['On effectue une phase de tests avec la fonction f1, qui atteint son minimum au point x=[1;1;1]']);

% Test n°1
disp(['*************************************************']);
x011 = [1;0;0]; 
disp(['Test pour f1, depuis le point de départ x011.']);
[x_min1,k1,flag1] = Newton(x011,@f1,@gradf1,@hessf1,1e-10,1e-10);
res = [res;{'f1',strjoin(string(x011)),strjoin(string(x_min1)),k1,strjoin(string(flag1))}];

% Test n°2
disp(['*************************************************']);
x012 = [10;3;-2.2]; 
disp(['Test pour f1, depuis le point de départ x012.']);
[x_min2,k2,flag2] = Newton(x012,@f1,@gradf1,@hessf1,1e-10,1e-10);
res = [res;{'f1',strjoin(string(x012)),strjoin(string(x_min2)),k2,strjoin(string(flag2))}];

% **** Phase 2 ****
disp(['*************************************************']);
disp(['On effectue une phase de tests avec la fonction f2, qui atteint son minimum au point x=[1;1]']);

% Test n°3
disp(['*************************************************']);
x021 = [-1.2;1]; 
disp(['Test pour f2, depuis le point de départ x021.']);
[x_min3,k3,flag3] = Newton(x021,@f2,@gradf2,@hessf2,1e-10,1e-10);
res = [res;{'f2',strjoin(string(x021)),strjoin(string(x_min3)),k3,strjoin(string(flag3))}];

% Test n°4
disp(['*************************************************']);
x022 = [10;0]; 
disp(['Pour f2, depuis le point de départ x022, on trouve : ']);
[x_min4,k4,flag4] = Newton(x022,@f2,@gradf2,@hessf2,1e-10,1e-10);
res = [res;{'f2',strjoin(string(x022)),strjoin(string(x_min4)),k4,strjoin(string(flag4))}];

% Test n°5
disp(['*************************************************']);
x023 = [0;1/200 + 1/(10^12)]; 
disp(['Test pour f2, depuis le point de départ x023.']);
[x_min5,k5,flag5] = Newton(x023,@f2,@gradf2,@hessf2,1e-10,1e-10);
res = [res;{'f2',strjoin(string(x023)),strjoin(string(x_min5)),k5,strjoin(string(flag5))}];



%% Construction du tableau des resultats 
cres = cell2table(res);
cres.Properties.VariableNames = {'Fonction','Point_de_depart','Point_minimisant','Nb_iterations','flag'};
writetable(cres,'testsNewton.xls','Sheet',1,'Range','A1');
disp(cres)



%% **** Definition des fonctions des test ****

% Definition de la fonction f1, de son gradient et sa hessienne 

function [res] = f1(x)

    x1 = x(1);
    x2 = x(2);
    x3 = x(3);
    
    res = 2*(x1+x2+x3-3)^2 + (x1-x2)^2 + (x2-x3)^2;

end

function [gf] = gradf1(f,x)

gf = 4*(x(1)+x(2)+x(3)-3)*ones(3,1) +[ 2*(x(1)-x(2)) ; -2*x(1)+4*x(2)-2*x(3) ; -2*(x(2)-x(3))];

end 

function [hf] = hessf1(f,x)

hf = [6 2 4 ; 2 8 2 ; 4 2 6];

end 

% Definition de la fonction f2, de son gradient et sa hessienne 

function [x] = f2(x1,x2)

x = 100*(x(2)-x(1)^2)^2 + (1-x1)^2;

end

function [gf] = gradf2(f,x)

gf = [ -400*(x(1)*x(2)-x(1)^3) - 2*(1-x(1)) ; 200*(x(2)-x(1)^2)];

end 

function [hf] = hessf2(f,x)

hf = [1200*x(1)^2 - 400*x(2) + 2 , -400*x(1) ; -400*x(1) , 200 ];

end 
