%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [lambda,numIter] = NonLinearNewton(lMin,lMax,phi,dPhi)
%       Résout l'équation phi(x) = 0 avec phi une fonction non lineaire de
%       la varibale x; phi : R -> R.

%% Parametres 
%       -lMin : borne inférieure pour la recherche dichotomique 
%       -lMax : borne supérieure pour la recherhe dichotomique
%       -phi : fonction non linéaire de la varibale x
%       -dPhi : fonction derivee de la fonction phi

%% Preconditions 
%       On note x* le point de R tel que phi(x*) = 0
%       1) lMin <= x*
%       2) x* <= lMax

%% Retours
%       -lambda : approximation d'un zero de la fonction f entre lMin et
%       lMax
%       -numIter : nombre d'itérations effectuées
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lambda,numIter] = NonLinearNewton(lMin,lMax,phi,dPhi)

numIter = 1;
epsilon = 10^(-4);


if (min(abs(phi(lMin)),abs(phi(lMax))) < epsilon)
    lambda = min(lMin,lMax); %permet de conserver phi(lMin)*phi(lMax) <= 0 et lambda petit petit
else 
    lambda = lMax;
end 

while ((min(abs(phi(lMin)),abs(phi(lMax))) >= epsilon) && (abs(phi(lambda)) >= epsilon))
    % Iteration de Newton
    derPhi = dPhi(lambda);
    lambdaN = lambda - (phi(lambda)/derPhi);
    if ((lMin <= lambdaN) && (lambdaN <= lMax) && abs((phi(lambdaN)) < 0.5*abs(phi(lambda))))
        % Cas où l'iteree est acceptee
        lambda = lambdaN;
        flag = 1;
    else 
        % Recherche d'un nouveau lambda par dichotomie $
        lambdaD = (lMin + lMax)/2;
        if (phi(lambdaD)*phi(lMax) <= 0)
            lMin = lambdaD;
        else 
            lMax = lambdaD;
        end 
        lambda = lambdaD; 
    end 
    
    numIter = numIter + 1;
end
        
    
