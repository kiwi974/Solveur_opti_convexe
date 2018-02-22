%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [xk,k,flag] = Newton(x0,f,gradf,hessf,tol1,tol2)
%       Cherche le point minimisant f à partir du point x0 grâce à des
%       approximations de type Newton 

%% Parametres 
%       -x0 : point de depart de la recherche 
%       -f : fonction que l'on cherche a minimiser
%       -gradf : fonction gradient de la fonction f
%       -hessf : fonction hessienne de la fonction f
%       -tol1 : precision pour la contrainte sur la stationnarite du
%       gradient
%       -tol2 : precision pour la contrainte sur la stationnarite du point
%       minimisant que l'on recherche 

%% Retours
%       -xk : point minimisant trouve, derniere valeur de ce point si l'on
%       a arrêté l'algortihme car on a effectue le nombre maximum
%       d'iterations
%       -k : nombre ditérations effectuées pour trouver xk
%       -flag : 1 -> stationnarite du gradient a la précision tol1
%               2 -> nombre d'iterations maximal atteint
%               3 -> stationnarite de xk d'une iteration a l'autre 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xk,k,flag] = Newton(x0,f,gradf,hessf,tol1,tol2)

%******************%
% Champ de donnees %
%******************%  

c = 10^(-2);
kmax = 10;
k = 1;

% Calcul du gradient au point de depart
xk_1 = x0;
gf0 = gradf(f,x0);
hf0 = hessf(f,x0);

% Calcul du premier point
dk = hf0\(-gf0);
xk = xk_1 + dk;
gf = gradf(f,xk);
hf = hessf(f,xk); 

% Tant qu'il n'y a pas convergence, on calcule un nouveau deplacement
while (norm(gf,2) > tol1*(norm(gf0,2) + c)) && (k <= kmax) && (norm(xk-xk_1,2) > tol2*(norm(xk,2) + c))
    xk = xk + (hf\(-gf));
    gf = gradf(f,xk);
    hf = hessf(f,xk);
    k = k + 1;
end 

if (norm(gf,2) <= tol1*(norm(gf0,2) + c))
    flag = 1
else 
    if (k > kmax)
        flag = 2
    else 
        flag = 3
    end
end


end 
