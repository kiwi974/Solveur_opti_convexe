%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [xk,k,flag] = RC(x0,f,gradf,hessf,delta,tol1,tol2,choix)
%       Calcul du minimum de la fonction f, à partir du point x0

%% Parametres 
%       x0 : point de depart de la recherche 
%       f : fonction que l'on cherche à minimiser (f est à valeur dans R)
%       graf : fonction gradient de la fonction f 
%       hessf : fonction hessienne de la fonction f
%       delta : rayon initial de la region de confiance
%       tol1 : tolerance pour la stationnarite de gradient
%       tol2 : tolerance pour la stationnarite du point x recherche
%       choix : 1 -> pas de Cauchy 
%               2 -> pas de More-Sorensen

%% Retour
%       xk : point minimisant f dans son domaine de definition
%       k : nombre d'iterations effectue
%       flag : 1 -> arrêt car stationnarite du gradient 
%              2 -> arrêt car nombre maximal d'itérations atteint
%              3 -> arrêt car stationnarite de l'iteree de x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [xk,k,flag] = RC(x0,f,gradf,hessf,delta,tol1,tol2,choix)

%choix = 1 -> Cauchy 
%choix = 2 -> More-Sorensen

%******************%
% Champ de donnees %
%******************%  
c = 10^(-10);
kmax = 20000;
k = 1;

eta1 = 0.1;
eta2 = 0.8;

gamma1 = 0.7;
gamma2 = 1.2;

% Calcul du gradient au point de depart
xk_1 = x0;
xk = x0;
norm_gf0 = norm(gradf(x0),2);
gk = gradf(xk);

% Initialisation de xk
xk_aux = Inf;

% Tant qu'il n'y a pas convergence, on calcule un nouveau deplacement
while (norm(gk,2) > tol1*(norm_gf0 + c)) & (k <= kmax) & (norm(xk_aux-xk_1,2) > tol2*norm(xk_1,2))
    % Calcul du modèle mk %
    fk = f(xk);
    gk = gradf(xk);
    hk = hessf(xk);
    % Resolution du sous-probleme approche
    if (choix == 1)
        sk = Cauchy(delta,gk,hk);
    else 
        [sk,~] = MoreSorensen(gk,hk,delta,1e-10);
    end 
    % Calcul des variaions % 
    p = (f(xk+sk)-fk)/((gk'*sk + (sk'*hk*sk)/2));
    % Mise a jour du point et de la region de confiance ù
    if (p >= eta2)
        delta = gamma2*delta;
        xk_1 = xk;
        xk = xk + sk;
        xk_aux = xk;
    else 
        if (p < eta1)
            delta = gamma1*delta;
        else 
            xk_1 = xk;
            xk = xk + sk;
            xk_aux = xk;
        end
    end      
    k = k + 1;   
end 




if (norm(gk,2) <= tol1*(norm_gf0 + c))
    %disp(['RC stoppée car (norm(gk,2) <= tol1*(norm(gf0,2) + c))']);
    flag = 1;
else 
    if (k > kmax)
        %disp(['RC stoppée car k > kmax']);
        flag = 2;
    else 
        %disp(['RC stoppée car (norm(xk-xk_1,2) < tol2)']);
        flag = 3;
    end 
end 
        


end 
