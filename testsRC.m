% Programme de test de l'algorithme de Newton avec régions de confiance 

res = {};

tol1 = 1e-8;
tol2 = 1e-8;
delta = 1;

% **** Phase 1 ****
disp(['*************************************************']);
disp(['On effectue une phase de tests avec la fonction f1, qui atteint son minimum au point x=[1;1;1]']);
disp(['*************************************************']);

% Test n°1
disp(['*************************************************']);
x011 = [1;0;0]; 
disp(['Test pour f1, depuis le point de départ x011.']);
t = tic;
[x_min1_C,k1_C,flag1_C] = RC(x011,@f1,@gradf1,@hessf1,delta,tol1,tol2,1);
tf = toc(t);
res = [res;{'f1',strjoin(string(x011)),delta,strjoin(string(x_min1_C)),k1_C,'Cauchy',flag1_C,tf}];
t = tic;
[x_min1_MS,k1_MS,flag1_MS] = RC(x011,@f1,@gradf1,@hessf1,delta,tol1,tol2,2);
tf = toc(t);
res = [res;{'f1',strjoin(string(x011)),delta,strjoin(string(x_min1_MS)),k1_MS,'MS',flag1_MS,tf}];


% Test n°2
disp(['*************************************************']);
x012 = [10;3;-2.2]; 
disp(['Test pour f1, depuis le point de départ x012.']);
t = tic;
[x_min2_C,k2_C,flag2_C] = RC(x012,@f1,@gradf1,@hessf1,delta,tol1,tol2,1);
tf = toc(t);
res = [res;{'f1',strjoin(string(x012)),delta,strjoin(string(x_min2_C)),k2_C,'Cauchy',flag2_C,tf}];
t = tic;
[x_min2_MS,k2_MS,flag2_MS] = RC(x012,@f1,@gradf1,@hessf1,delta,tol1,tol2,2);
tf = toc(t);
res = [res;{'f1',strjoin(string(x012)),delta,strjoin(string(x_min2_MS)),k2_MS,'MS',flag2_MS,tf}];

% **** Phase 2 ****
disp(['*************************************************']);
disp(['On effectue une phase de tests avec la fonction f2, qui atteint son minimum au point x=[1;1]']);
disp(['*************************************************']);

% Test n°3
disp(['*************************************************']);
x021 = [-1.2;1]; 
disp(['Test pour f2, depuis le point de départ x021.']);
t = tic;
[x_min3_C,k3_C,flag3_C] = RC(x021,@f2,@gradf2,@hessf2,delta,tol1,tol2,1);
tf = toc(t);
res = [res;{'f2',strjoin(string(x021)),delta,strjoin(string(x_min3_C)),k3_C,'Cauchy',flag3_C,tf}];
t = tic;
[x_min3_MS,k3_MS,flag3_MS] = RC(x021,@f2,@gradf2,@hessf2,delta,tol1,tol2,2);
tf = toc(t);
res = [res;{'f2',strjoin(string(x021)),delta,strjoin(string(x_min3_MS)),k3_MS,'MS',flag3_MS,tf}];

% Test n°4
disp(['*************************************************']);
x022 = [10;0]; 
disp(['Test pour f2, depuis le point de départ x022.']);
t = tic;
[x_min4_C,k4_C,flag4_C] = RC(x022,@f2,@gradf2,@hessf2,delta,tol1,tol2,1);
tf = toc(t);
res = [res;{'f2',strjoin(string(x022)),delta,strjoin(string(x_min4_C)),k4_C,'Cauchy',flag4_C,tf}];
t = tic;
[x_min4_MS,k4_MS,flag4_MS] = RC(x022,@f2,@gradf2,@hessf2,delta,tol1,tol2,2);
tf = toc(t);
res = [res;{'f2',strjoin(string(x022)),delta,strjoin(string(x_min4_MS)),k4_MS,'MS',flag4_MS,tf}];

% Test n°5
disp(['*************************************************']);
x023 = [0;1/200 + 1/(10^12)]; 
disp(['Test pour f2, depuis le point de départ x023.']);
t = tic;
[x_min5_C,k5_C,flag5_C] = RC(x023,@f2,@gradf2,@hessf2,delta,tol1,tol2,1);
tf = toc(t);
res = [res;{'f2',strjoin(string(x023)),delta,strjoin(string(x_min5_C)),k5_C,'Cauchy',flag5_C,tf}];
t = tic;
[x_min5_MS,k5_MS,flag5_MS] = RC(x023,@f2,@gradf2,@hessf2,delta,tol1,tol2,2);
tf = toc(t);
res = [res;{'f2',strjoin(string(x023)),delta,strjoin(string(x_min5_MS)),k5_MS,'MS',flag5_MS,tf}];


%% Construction du tableau des resultats 

cres = cell2table(res);
cres.Properties.VariableNames = {'Fonction','Point_depart','delta','Point_minimiant','Nb_iter','Methode_dans_RC','flag','Temps_exec'};
writetable(cres,'testsRC.xls','Sheet',1,'Range','A1');
disp(cres)


%% Tests pour les temps d'execution
% On prend volontairement un point éloigné du minimum de chaque fonction
% pour voir quelle méthode l'atteint le plus vite

resTime = {};

tol1 = 1e-12;
tol2 = 1e-12;

% Test n°6
disp(['*************************************************']);
x100 = [100;100;100]; 
disp(['Test pour f1, depuis le point de départ x100.']);
t = tic;
[x_min6_C,k6_C,flag6_C] = RC(x100,@f1,@gradf1,@hessf1,delta,tol1,tol2,1);
tf = toc(t);
resTime = [resTime;{'f1',strjoin(string(x100)),delta,strjoin(string(x_min6_C)),k6_C,'Cauchy',flag6_C,tf}];
t = tic;
[x_min6_MS,k6_MS,flag6_MS] = RC(x100,@f1,@gradf1,@hessf1,delta,tol1,tol2,2);
tf = toc(t);
resTime = [resTime;{'f1',strjoin(string(x100)),delta,strjoin(string(x_min6_MS)),k6_MS,'MS',flag6_MS,tf}];

% Test n°7
disp(['*************************************************']);
x100_2 = [100;100]; 
disp(['Test pour f2, depuis le point de départ x100_2.']);
t = tic;
[x_min7_C,k7_C,flag7_C] = RC(x100_2,@f2,@gradf2,@hessf2,delta,tol1,tol2,1);
tf = toc(t);
resTime = [resTime;{'f2',strjoin(string(x100_2)),delta,strjoin(string(x_min7_C)),k7_C,'Cauchy',flag7_C,tf}];
t = tic;
[x_min7_MS,k7_MS,flag7_MS] = RC(x100_2,@f2,@gradf2,@hessf2,delta,tol1,tol2,2);
tf = toc(t);
resTime = [resTime;{'f2',strjoin(string(x100_2)),delta,strjoin(string(x_min7_MS)),k7_MS,'MS',flag7_MS,tf}];



%% Construction du tableau des resultats sur les temps d'execution

cresTime = cell2table(resTime);
cresTime.Properties.VariableNames = {'Fonction','Point_depart','delta','Point_minimiant','Nb_iter','Methode_dans_RC','flag','Temps_exec'};
writetable(cresTime,'testsRC_execution_time.xls','Sheet',1,'Range','A1');
disp(cresTime)

