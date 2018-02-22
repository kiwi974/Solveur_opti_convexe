% Programme de test de l'algorithme de Moré-Sorensen

epsZero = 10^(-10);

res = {};

% Test 1
disp(['*****************************************************************']);

disp(['Test n°1 : quadratique 1']);
g1 = [0 ; 0];
H1 = [7 0 ; 0 2];
[s1,lambda1,flag1] = MoreSorensen(g1,H1,0.7,epsZero);
res = [res;{'q1',0.7,strjoin(string(s1)),lambda1,flag1}];

% Test 2
disp(['*****************************************************************']);

disp(['Test n°2 : quadratique 2']);
g2 = [6 ; 2];
H2 = [7 0 ; 0 2];
[s2,lambda2,flag2] = MoreSorensen(g2,H2,0.7,epsZero);
res = [res;{'q2',0.7,strjoin(string(s2)),lambda2,flag2}];

% Test 3
disp(['*****************************************************************']);

disp(['Test n°3 : quadratique 3']);
g3 = [-2 ; 1];
H3 = [-2 0 ; 0 10];
[s3,lambda3,flag3] = MoreSorensen(g3,H3,0.7,epsZero);
res = [res;{'q3',0.7,strjoin(string(s3)),lambda3,flag3}];

% Test 4
disp(['*****************************************************************']);

disp(['Test n°4 : quadratique 4']);
g4 = [0 ; 0];
H4 = [-2 0 ; 0 10];
[s4,lambda4,flag4] = MoreSorensen(g4,H4,1,epsZero);
res = [res;{'q4',1,strjoin(string(s4)),lambda4,flag4}];

% Test 5
disp(['*****************************************************************']);

disp(['Test n°5 : quadratique 5']);
g5 = [2 ; 3];
H5 = [4 6 ; 6 5];
[s5,lambda5,flag5] = MoreSorensen(g5,H5,1,epsZero);
res = [res;{'q5',1,strjoin(string(s5)),lambda5,flag5}];

% Test 6
disp(['*****************************************************************']);

disp(['Test n°6 : quadratique 6']);
g6 = [2 ; 0];
H6 = [4 0 ; 0 -15];
[s6,lambda6,flag6] = MoreSorensen(g6,H6,1,epsZero);
res = [res;{'q6',1,strjoin(string(s6)),lambda6,flag6}];


%% Construction du tableau des resultats 
cres = cell2table(res);
cres.Properties.VariableNames = {'Quadratique','Rayon_RC','Pas_de_MS','Multiplicateur','flag'};
writetable(cres,'testsMoreSorensen.xls','Sheet',1,'Range','A1');
disp(cres)