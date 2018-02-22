%% Script de test de la fonction effectuant la methode du Lagrangien
% augmenté


% Definition des constantes de depart conformément au sujet 

alpha = 0.1;
beta = 0.9;
eta = 0.1;
eta_chap = 0.1258925;
mu = (eta_chap/eta)^(1/alpha);
tau = 1.4;
lambda = 0;
delta = 1;

% Definition des points de depart
x0_11 = [0;1;1];
x0_12 = [0.5;1.25;1];
x0_21 = [1;0];
u = sqrt(3)/2;
x0_22 = [u;u];

%% Test de l'influence du choix de la méthode 



%% Tests de l'influence du choix de lambda au départ avec f1

% Etude de l'influence du parametre lambda avec alpha = 0.1, eta = 0.1,
% eta_chap = 0.1258925 (donc mu fixe), et le pas de Cauchy

res_lambda = {};
lambda_val = [-6 -2 0 0.1 0.5 1 1.2 1.8 2 3 4 8 12 20 30 40 60];
for i =1:length(lambda_val)
    disp(['Essai de la valeur ' num2str(lambda_val(i)) ' de lambda'])
    [x_res,l_res,mu_res,iter_res,f_res] = LA(x0_11,lambda_val(i),mu,alpha,beta,eta_chap,tau,delta,@LA1,@gradLA1,@hessLA1,@c1,1);
    res_lambda = [res_lambda;{'f1',strjoin(string(x0_11)),lambda_val(i),strjoin(string(x_res)),l_res,mu_res,iter_res,f_res}];
end 
cres_lambda = cell2table(res_lambda);
cres_lambda.Properties.VariableNames = {'Fonction','Point_de_depart','lambda_0','Point_minimisant','lambda_res','mu_res','Nb_iterations','flag'};
writetable(cres_lambda,'testsLA_lambda.xls','Sheet',1,'Range','A1');
disp(cres_lambda)






%% Tests de l'influence du choix de lambda au départ avec f2

% Etude de l'influence du parametre lambda avec alpha = 0.1, eta = 0.1,
% eta_chap = 0.1258925 (donc mu fixe), et le pas de Cauchy

res_lambda2 = {};
lambda_val = [-6 -2 0 0.1 0.5 1 1.2 1.8 2 3 4 8 12 20 30 40 60];
for i =1:length(lambda_val)
    disp(['Essai de la valeur ' num2str(lambda_val(i)) ' de lambda'])
    [x_res,l_res,mu_res,iter_res,f_res] = LA(x0_21,lambda_val(i),mu,alpha,beta,eta_chap,tau,delta,@LA2,@gradLA2,@hessLA2,@c2,1);
    res_lambda2 = [res_lambda2;{'f2',strjoin(string(x0_21)),lambda_val(i),strjoin(string(x_res)),l_res,mu_res,iter_res,f_res}];
end 
cres_lambda2 = cell2table(res_lambda2);
cres_lambda2.Properties.VariableNames = {'Fonction','Point_de_depart','lambda_0','Point_minimisant','lambda_res','mu_res','Nb_iterations','flag'};
writetable(cres_lambda2,'testsLA_lambda2.xls','Sheet',1,'Range','A1');
disp(cres_lambda2)





%% Tests de l'influence du choix de tau au départ avec f1

% Etude de l'influence du parametre tau avec alpha = 0.1, eta = 0.1,
% eta_chap = 0.1258925 (donc mu fixe), et lambda = 1

res_tau = {};
tau_val = [1.01 1.5 1.7 2 2.5 3 6 10 15 20 30];
for i =1:length(tau_val)
    disp(['Essai de la valeur ' num2str(tau_val(i)) ' de tau'])
    [x_res,l_res,mu_res,iter_res,f_res] = LA(x0_11,lambda,mu,alpha,beta,eta_chap,tau_val(i),delta,@LA1,@gradLA1,@hessLA1,@c1,1);
    res_tau = [res_tau;{'f1',strjoin(string(x0_11)),tau_val(i),strjoin(string(x_res)),l_res,mu_res,iter_res,f_res}];
end 
cres_tau = cell2table(res_tau);
cres_tau.Properties.VariableNames = {'Fonction','Point_de_depart','tau_0','Point_minimisant','lambda_res','mu_res','Nb_iterations','flag'};
writetable(cres_tau,'testsLA_tau.xls','Sheet',1,'Range','A1');
disp(cres_tau)


%% Tests de l'influence du choix de tau au départ avec f2

% Etude de l'influence du parametre tau avec alpha = 0.1, eta = 0.1,
% eta_chap = 0.1258925 (donc mu fixe), et lambda = 1

res_tau2 = {};
tau_val = [1.01 1.5 1.7 2 2.5 3 6 10 15 20 30];
for i =1:length(tau_val)
    disp(['Essai de la valeur ' num2str(tau_val(i)) ' de tau'])
    [x_res,l_res,mu_res,iter_res,f_res] = LA(x0_21,lambda,mu,alpha,beta,eta_chap,tau_val(i),delta,@LA2,@gradLA2,@hessLA2,@c2,1);
    res_tau2 = [res_tau2;{'f2',strjoin(string(x0_21)),tau_val(i),strjoin(string(x_res)),l_res,mu_res,iter_res,f_res}];
end 
cres_tau2 = cell2table(res_tau2);
cres_tau2.Properties.VariableNames = {'Fonction','Point_de_depart','tau_0','Point_minimisant','lambda_res','mu_res','Nb_iterations','flag'};
writetable(cres_tau2,'testsLA_tau2.xls','Sheet',1,'Range','A1');
disp(cres_tau2)







%% Definition du lagrangien et des derivees de chacun des deux tests 

% LA1

function [c] = c1(x)
    c = x(1) + x(3) - 1;
end

function [y] = LA1(x,lambda,mu)
    y = f1(x) + lambda*(x(1)+x(3)-1) +mu*(x(1)+x(3)-1)^2;
end

function [y] = gradLA1(x,lambda,mu)
    y = gradf1(x) + (lambda+(mu/2)*(x(1)+x(3)-1))*[1;0;1];
end

function [y] = hessLA1(x,lambda,mu)
    y = [7 2 5 ; 2 8 2 ; 5 2 7];
end
    
% -----------LA2---------------

% Contrainte

function [c] = c2(x)
    c = x(1)^2 + x(2)^2 - 1.5;
end

function [gc] = gradc2(x)
    gc = [2*x(1) ; 2*x(2)];
end

function [hc] = hessc2(x)
    hc = [2 0 ; 0 2];
end

% Lagrangien

function [y] = LA2(x,lambda,mu)
    c = c2(x);
    y = f2(x) + lambda*c +(mu/2)*c^2;
end

function [y] = gradLA2(x,lambda,mu)
    c = c2(x);
    gc = gradc2(x);
    y = gradf2(x) + lambda*gc + mu*gc*c;
end

function [y] = hessLA2(x,lambda,mu)
    c = c2(x);
    gc = gradc2(x);
    hc = hessc2(x);
    y = hessf2(x) + (lambda + mu*c)*hc + gc*gc';
end
