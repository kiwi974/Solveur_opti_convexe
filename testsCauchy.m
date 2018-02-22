% Tests de la fonction trouvant le pas de Cauchy 

% Test n°1

disp(['*************************************************************']);
disp(['Test n°1']);
disp(['*************************************************************']);
disp(['Avec la quadratique 1 et un rayon delta=1, on trouve le pas de Cauchy : ']);
g1 = [0;0];
H1 = [7 0 ; 0 2];

s1 = Cauchy(1,g1,H1)

% Test n°2

disp(['*************************************************************']);
disp(['Test n°2']);
disp(['*************************************************************']);
disp(['Avec la quadratique 2 et un rayon delta=1, on trouve le pas de Cauchy : ']);
g2 = [6;2];
H2 = [7 0 ; 0 2];

s2 = Cauchy(1,g2,H2)

% Test n°3

disp(['*************************************************************']);
disp(['Test n°3']);
disp(['*************************************************************']);
disp(['Avec la quadratique 2 et un rayon delta=1, on trouve le pas de Cauchy : ']);
g2 = [6;2];
H2 = [7 0 ; 0 2];
s3 = Cauchy(10^(-2),g2,H2)

% Test n°4

disp(['*************************************************************']);
disp(['Test n°3']);
disp(['*************************************************************']);
disp(['Avec la quadratique 3 et un rayon delta=1, on trouve le pas de Cauchy : ']);
g3 = [-2 ; 1];
H3 = [-2 0 ; 0 10];
s4 = Cauchy(1,g3,H3)

% Test n°5

disp(['*************************************************************']);
disp(['Test n°5']);
disp(['*************************************************************']);
disp(['Avec la quadratique 3 et un rayon delta=10, on trouve le pas de Cauchy : ']);
g3 = [-2 ; 1];
H3 = [-2 0 ; 0 10];
s5 = Cauchy(10,g3,H3)

%% Definition des quadratiques de test 

% Quadratique n°1
g1 = [0;0];
H1 = [7 0 ; 0 2];

% Quadratique n°2
g2 = [6;2];
H2 = [7 0 ; 0 2];

% Quadratique 3
g3 = [-2 ; 1];
H3 = [-2 0 ; 0 10];

