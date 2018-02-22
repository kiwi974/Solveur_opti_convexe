% Programme de test de l'algorithme de la recherche de 0 par la méthode de
% Newton pour des fonctions d'une variable réelle à valeurs réelles 

res = {};

% Test 1
disp(['*****************************************************************']);

disp(['Test n°1 : a = 4 ; b = 2 ; c = 36 ; d = 14 ; delta = 0.5']);
t = tic;
[lambda1,n1] = NonLinearNewton(-8,10,@phi1,@dphi1);
tfin = toc(t);
res = [res;{'phi1',-8,10,lambda1,n1,phi1(lambda1)}];

disp(['*****************************************************************']);

% Test 2 
disp(['Test n°2 : a = 4 ; b = -38 ; c = 400 ; d = 20 ; delta = 0.2']);
t = tic;
[lambda2,n2] = NonLinearNewton(80,100,@phi2,@dphi2);
tfin = toc(t);
res = [res;{'phi2',80,100,lambda2,n2,phi2(lambda2)}];

disp(['*****************************************************************']);

% Test 3

disp(['Test n°3 : a = 4 ; b = -38 ; c = 400 ; d = 20 ; delta = 0.7']);
t = tic;
[lambda3,n3] = NonLinearNewton(30,55,@phi3,@dphi3);
tfin = toc(t);
res = [res;{'phi3',30,55,lambda3,n3,phi3(lambda3)}];

disp(['*****************************************************************']);


% Test 4
% On teste en cherchant en plus les bornes de l'intervalle de recherche 

disp(['Test n°4 : On teste en cherchant en plus les bornes de l intervalle de recherche ']);
t = tic;
[lmin,lmax] = bornesDicho(@phi3,37);
[lambda4,n4] = NonLinearNewton(lmin,lmax,@phi3,@dphi3);
tfin = toc(t);
res = [res;{'phi3',lmin,lmax,lambda4,n4,phi3(lambda4)}];

disp(['*****************************************************************']);


%% Contruction du tableau des resultats 

cres = cell2table(res);
cres.Properties.VariableNames = {'Fonction','lMin','lMax','Zero','Nb_iter','Phi_zero_approxime'};
writetable(cres,'testsNewton.xls','Sheet',1,'Range','A1');
disp(cres)



%% Fonctions de test 

% Babaris fonction 

function [x] = phiN(lambda,a,b,c,d,delta)

    x = a/((lambda+b)^2) + c/((lambda+d)^2) - delta^2;
 
end

function [x] = phiNINV(lambda,a,b,c,d,delta)

    x = 1/(a/((lambda+b)^2) + c/((lambda+d)^2)) - 1/(delta^2);
 
end

function [x] = d_phiN(lambda,a,b,c,d)

    x = (-2*a)/((lambda+b)^3) + (-2*c)/((lambda+d)^3);
 
end

function [x] = d_phiNINV(lambda,a,b,c,d)

    x = ((2*a)/((lambda+b)^3) + (2*c)/((lambda+d)^3))/((a/((lambda+b)^2) + c/((lambda+d)^2))^2);
 
end

%**************************************************************************
% Fonction test 1
function [z] = phi1(lambda)
    z = phiNINV(lambda,4,2,36,14,0.5);
end 

function [z] = dphi1(lambda)
    z = d_phiNINV(lambda,4,2,36,14);
end
%**************************************************************************

%**************************************************************************
% Fonction test 2
function [z] = phi2(lambda)
    z = phiNINV(lambda,4,-38,400,20,0.2);
end 


function [z] = dphi2(lambda)
    z = d_phiNINV(lambda,4,-38,400,20);
end
%**************************************************************************

%**************************************************************************
% Fonction test 3
function [z] = phi3(lambda)
    z = phiNINV(lambda,4,-38,400,20,0.7);
end 

function [z] = dphi3(lambda)
    z = d_phiNINV(lambda,4,-38,400,20);
end
%**************************************************************************



% Fonction pour trouver lmin et lmax pour utiliser la dichotomie dans
% NonLinearNewtion

function [lmin,lmax] = bornesDicho(f,lambda1)
    %lmin = max(0,-lambda1) + 10^(-8);
    lmin = 29;
    while (f(lmin) > 0) 
        lmin = lmin + 0.5;
    end
    lmax = lmin+0.01;
    while (f(lmax) < 0)
        lmax = lmax + 0.5;
    end
end
