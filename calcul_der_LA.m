%% Definition de la fonction LA et de ses derivees 

function [y] = LA1(x,lambda,mu,f,c)
	constraint = c(x);
	y = f(x) + lambda'*constraint + 0.5*mu*norm(constraint,2)^2;
end


% param : 
%	- gradc : tableau des gradients des fonctions composantes de c 

function [y] = gradLA1(x,lambda,mu,gradf,c,gradc)
	% Construction de la jacobienne de c à partir de gradc(x)
	jacoc = Jc(x);
	y = gradf(x) + lambda'*jacoc + mu*jacoc'*c(x);
end


% param : 
%	- gradc : tableau des gradients des fonctions composantes de c 
%	- hessc : tableau des hessiennes des fonctions composantes de c

function [y] = hessLA1(x,lambda,mu,hessf,c,gradc,hessci)
	% Evaluation des differentes fonctions au point x
	vectc = c(x);
	nablac = gradc(x);
	nablac2 = hessc(x);
	% Construction de la jacobienne de c à partir de gradc(x)
	% ???
	% Calcul des termes de la hessienne de LA
	sum1 = 0;
	sum2 = 0;
	for i=1:m
		sum1 = sum1 + lambda(i)*nablac(i);
		sum2 = sum2 +  vectc(i)*nablac2(i);	
	end
        % Calcul du resultat LA(x)
	y = hessf + sum1 + sum2 + jacoc'*jacoc;
end    