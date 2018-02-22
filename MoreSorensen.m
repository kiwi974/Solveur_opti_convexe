%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [s,lambda,flag] = MoreSorensen(g,H,delta,epsZero)
%       Fonction effectuant le pas de More-Sorensen pour une quadratique 
%       de terme de premier ordre g et de terme de second ordre H, avec un
%       rayon de la région de confiance delta

%       On note l1 la plus grand valeur propre de H

%% Parametres 
%       -g : terme de premier ordre de la quadratique
%       -H : terme de second ordre de la quadratique
%       -delta : rayon de la region de confiance
%       -epsZero : seuil pour verifier si le vecteur propre associe a l1 
%       est orthogonal à g

%% Retours
%       -s : pas de More-Sorensen
%       -lambda : multiplicateur associé dans le problème KKT
%       -flag : 0 -> Solution intérieure (H > 0 et pas de Cauchy)
%               1 -> Solution sur la frontière et le multiplicateur de
%               Lagrange associe n'est pas l1
%               2 -> Solution sur la frontière et le multiplicateur de
%               Lagrange associe est l1 : s à l'extérieur de la RC
%               3 -> Solution sur la frontière et le multiplicateur de
%               Lagrange associe n'est pas l1 : 'cas difficie'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [s,lambda,flag] = MoreSorensen(g,H,delta,epsZero)

%epsZero = 10^(-10);

% Dimenson de H (h est carrée)
n = length(H(1,:));

% Calcul du pas de Newton 
sN = H\(-g);

% Decomposition spectrale de H
[Q,D] = eig(H,'vector');  % H = Q*D*Q'
[D,index] = sort(D); % valeurs propres ordonnées par ordre croissant

Q = Q(:,index);

% Boolean valant vrai ssi H est définie positive
defPosH = isempty(find(D < 0));

% Cas d'une solution intérieure
if (defPosH && norm(sN,2)^2 <= delta^2)
    %disp(['H est semie-définie positive et le pas de Cauchy permet de conclure']);
    aux = find((D >0));
    l = aux(1);
    s = sN;
    lambda = 0;
    flag = 0;
    
else % Cas d'une solution sur la frontière 
    if (abs(Q(:,1)'*g) > epsZero) % Sous-cas où le multiplicateur de Lagrange associé n'est pas D(1)
        %disp(['Le multiplicateur de Lagrange associé n est pas D(1)']);
        [lmin,lmax] = bornesDicho(@(lambda)fun_s_norm_Newton(lambda,delta,n,Q,D,g,1),D(1));
        [lambda,~] = NonLinearNewton(lmin,lmax,@(lambda)fun_s_norm_Newton(lambda,delta,n,Q,D,g,1),@(lambda)fun_d_s_norm(lambda,delta,n,Q,D,g,1));
        s = fun_s(lambda,delta,n,Q,D,g,1);
        flag = 1;
    else % Sous-cas où le multiplicateur de Lagrange associé est D(1)
        %disp(['Le multiplicateur de Lagrange associé est D(1)']);
        sum_lambda1_norm = fun_s_norm(-D(1),delta,n,Q,D,g,2);
        if (norm(sum_lambda1_norm,2) > delta^2)
            %disp(['Cas ou la norme de s dépasse delta']);
            [lmin,lmax] = bornesDicho(@(lambda)fun_s_norm_Newton(lambda,delta,n,Q,D,g,2),D(1));
            %disp(['lmin,lmax OK ']);
            [lambda,~] = NonLinearNewton(lmin,lmax,@(lambda)fun_s_norm_Newton(lambda,delta,n,Q,D,g,2),@(lambda)fun_d_s_norm(lambda,delta,n,Q,D,g,2));
            %disp(['zero par Newton ok ']);
            s = fun_s(lambda,delta,n,Q,D,g,2);
            flag = 2;
        else % "Cas difficile"
            %disp(['Cas difficile']);
            % Pythagore pour trouver alpha
            lambda = -D(1);
            alpha = sqrt((delta^2 - sum_lambda1_norm)/(norm(Q(:,1),2)^2));
            s = fun_s(lambda,delta,n,Q,D,g,2) + alpha*Q(:,1);
            flag = 3;
        end
    end 
end

end



%% Définition des fonctions auxiliaires

% La fonction s et la fonction en norme 2 associée
function [s] = fun_s(lambda,delta,n,Q,D,g,indice)
    somme = 0;
    for i =indice:n
        somme = somme + ((Q(:,i)'*g)/(D(i)+lambda))*Q(:,i);
    end
    s = -somme;
end

function [s] = fun_s_norm(lambda,delta,n,Q,D,g,indice)
    somme = 0;
    for i =indice:n
        somme = somme + (Q(:,i)'*g)^2/((D(i)+lambda)^2);
    end
    s = somme;
end

% Dérivée de s et de sa fonction en norme 2 associée (pour trouver les
% zéros avec Newton)

function [s] = fun_d_s(lambda,delta,n,Q,D,g,indice)
    somme = 0;
    for i =indice:n
        somme = somme + 2*((Q(:,i)'*g)/((D(i)+lambda)^2))*Q(:,i);
    end
    s = somme;
end

function [s] = fun_d_s_norm(lambda,delta,n,Q,D,g,indice)
    somme = 0;
    for i =indice:n
        somme = somme + (Q(:,i)'*g)^2/((D(i)+lambda)^3);
    end
    s = -2*somme;
end

% Fonctions dont on veut trouver les zéros pour Newton 

function [s] = fun_s_norm_Newton(lambda,delta,n,Q,D,g,indice)
    somme = 0;
    for i =indice:n
        somme = somme + (Q(:,i)'*g)^2/((D(i)+lambda)^2);
    end
    s = somme - delta^2;
end


function [s] = fun_s_Newton(lambda,delta,n,Q,D,g,indice)
    s = 1/(fun_s_norm(lambda,delta,n,Q,D,g,indice)) - 1/(delta^2);
end

function [s] = fun_d_s_Newton(lambda,delta,n,Q,D,g,indice)
    sum = 0;
    for i =indice:n
        sum = sum + (Q(:,i)'*g)^2/((D(i)+lambda)^3);
    end
    s = 2*sum/((fun_s_norm(lambda,delta,n,Q,D,g,indice))^2);
end

% Fonction pour trouver lmin et lmax pour utiliser la dichotomie dans
% NonLinearNewtion

function [lmin,lmax] = bornesDicho(f,lambda1)
    lmin = max(0,-lambda1);
    %disp(['Recherche de lmin (depuis ' num2str(lmin) ') et lmax']);
    %disp(['f(lmin) = ' num2str(f(lmin))]);
    while (f(lmin) < 0) 
        lmin = lmin + 0.1;
        %disp(['lmin vaut ' num2str(lmin) ' et f(lmin) vaut ' num2str(f(lmin))]);
        %pause
    end
    lmax = lmin+0.01;
    while (f(lmax) > 0)
        lmax = lmax + 50;
        %disp(['lmax vaut ' num2str(lmax) ' et f(lmax) vaut ' num2str(f(lmax))]);
        %pause
    end
end
