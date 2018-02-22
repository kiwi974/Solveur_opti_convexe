%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [x,lambda,mu,nbIter,flag] = LA(x,lambda,mu,alpha,beta,eta_chap,tau,delta,fLA,gradLA,hessLA,c,choix)
%       Calcul du minimum de la fonction f, grâce à la méthode du
%       multiplicateur de Lagrange

%% Parametres 
%       x : point de depart de la recherche 
%       mu : pénalité initiale sur la norme des contraintes
%       alpha : permet d'initialiser eta (cf. sujet)
%       beta : permet de mettre à jour eta (cf. sujet)
%       eta_chap : permet le calcul de eta (cf.sujet)
%       tau : permet l'agrandissement du paramètre de pénalité
%       delta : rayon initial de la méthode des régions de confiance utilisée
%       fLA : lagrangien associé à f
%       grafLA : dérivée du lagrangien associé à f
%       hessf : dérivée seconde du lagrangien associé à f
%       delta : rayon initial de la region de confiance
%       c : fonction associée aux contraintes 
%       choix : 1 -> Cauchy dans la méthode des régions de confiance
%               2 -> More-Sorensen dans la méthode des régions de confiance

%% Retour
%       x : point minimisant f 
%       lambda : multiplicateur de Lagrange
%       mu : pénalité à la fin de l'algorithme
%       nbIter : nombre d'itérations réalisées 
%       flag : 1 -> arrêt car stationnarite du gradient et contraintes
%                   respectées
%              2 -> arrêt car nombre maximal d'itérations atteint
%              3 -> arrêt car stationnarite de l'iteree de x et contraintes
%                   respectées
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,lambda,mu,nbIter,flag] = LA(x,lambda,mu,alpha,beta,eta_chap,tau,delta,fLA,gradLA,hessLA,c,choix)


% Calcul des parametres calculables %
prec = 10^(-6);
epsilon = 1/mu;
eta = eta_chap/(mu^alpha);

% Nombre maximum d'itérations %
nbIterMax = 100;
nbIter = 1;

x_aux = Inf;
x_prec = x;
tol2 = 1e-8;
tol3 = 1e-8;

while (((norm(gradLA(x,lambda,0),2) > prec) || (norm(c(x),2)^2 > tol2)) && (nbIter <= nbIterMax) && ((norm(x_aux-x_prec,2) > tol3*norm(x_prec,2)) || (norm(c(x),2)^2 > tol2)))

    %disp(['Iteration n° ' num2str(nbIter)]);
    
    % Resolution approximative du problème sans contrainte %
    % Utilisation de l'algortihme des régions de confiance %
    x_prec = x;
    [x,~,~] = RC(x,@(xvar)fLA(xvar,lambda,mu),@(xvar)gradLA(xvar,lambda,mu),@(xvar)hessLA(xvar,lambda,mu),delta,epsilon,tol2,choix);
    x_aux = x;
    
    
    if (norm(c(x),2)^2 <= eta)
        lambda = lambda + mu*c(x);
        % -- mu est inchange
        epsilon = epsilon/mu;
        eta = eta/(mu^beta);
    else 
        % -- lambda est inchange
        mu = tau*mu;
        epsilon = epsilon/mu;
        eta = eta_chap/(mu^alpha);
    end 
    nbIter = nbIter + 1;	
    %disp(['Iteration : ' num2str(nbIter)])
end 



if ((norm(gradLA(x,lambda,0),2) <= prec) && (norm(c(x),2)^2 <= tol2))
    flag = 1;
else 
    if (nbIter > nbIterMax)
        flag = 2;
    else 
        flag = 3;
    end 
end

end 




    
   



  
