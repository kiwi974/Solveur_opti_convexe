%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% s = Cauchy(deltak,g,H)
%       Calcul du pas de Cauchy pour une fonction f de terme de premier ordre
%       g, et terme de second ordre H

%% Parametres 
%       -deltak : rayon de la region de recherche
%       -g : terme de premier ordre de la fonction f
%       -H : terme de second ordre de la fonction f

%% Retours
%       -s : pas de Cauchy trouve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = Cauchy(deltak,g,H)

%Calcul du terme du second ordre%

norm_grad = norm(g,2);

if (norm_grad ~= 0)

	curvature = g'*H*g;

	tsd = (norm_grad^2)/curvature;

	limit = deltak/norm_grad;

	t = 0;

	if (curvature > 0)
		if (tsd < limit) 
			t = tsd;
		else
			t = limit;
		end
	else    
		t = limit;
	end

	% t permet alors de trouver la solution au problème de Cauchy%
	s = -t*g;

    %if (curvature > 0)
	%	disp(['curvature > 0 et limit = ' num2str(limit) ' et tsd = ' num2str(tsd), ' donc t = ' num2str(t)]);
    %else  
	%	disp(['curvature <= 0 et limit = ' num2str(limit) ' donc t = ' num2str(t) ]);
    %end




else 
	%disp(['!!Attention, le gradient est nul!!']);
	s = 0*g;
end

        
        
    
    
    
