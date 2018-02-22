res = {};

% Test n°1
disp(['******************************************']);
disp(['Test n°1 : f1 et point réalisable']);

[x1,l1,m1,n1,f1] = LA(x0_11,lambda,mu,alpha,beta,eta_chap,tau,delta,@LA1,@gradLA1,@hessLA1,@c1,1)
res = [res;{'LA1',strjoin(string(x0_11)),strjoin(string(x1)),11,m1,n1,'Cauchy',f1}];

[x2,l2,m2,n2,f2] = LA(x0_11,lambda,mu,alpha,beta,eta_chap,tau,delta,@LA1,@gradLA1,@hessLA1,@c1,2)
res = [res;{'LA1',strjoin(string(x0_11)),strjoin(string(x2)),12,m2,n2,'MS',f2}];

disp(['******************************************']);

% Test n°1
disp(['******************************************']);
disp(['Test n°2 : f1 et point non réalisable']);

[x3,l3,m3,n3,f3] = LA(x0_12,lambda,mu,alpha,beta,eta_chap,tau,delta,@LA1,@gradLA1,@hessLA1,@c1,1)
res = [res;{'LA1',strjoin(string(x0_12)),strjoin(string(x3)),13,m3,n3,'Cauchy',f3}];

[x4,l4,m4,n4,f4] = LA(x0_12,lambda,mu,alpha,beta,eta_chap,tau,delta,@LA1,@gradLA1,@hessLA1,@c1,2)
res = [res;{'LA1',strjoin(string(x0_12)),strjoin(string(x4)),14,m4,n4,'MS',f4}];

disp(['******************************************']);

% Test n°3
disp(['******************************************']);
disp(['Test n°3 : f2 et point non réalisable']);

[x5,l5,m5,n5,f5] = LA(x0_21,lambda,mu,alpha,beta,eta_chap,tau,delta,@LA2,@gradLA2,@hessLA2,@c2,1)
res = [res;{'LA2',strjoin(string(x0_21)),strjoin(string(x5)),15,m5,n5,'Cauchy',f5}];

[x6,l6,m6,n6,f6] = LA(x0_21,lambda,mu,alpha,beta,eta_chap,tau,delta,@LA2,@gradLA2,@hessLA2,@c2,2);
res = [res;{'LA2',strjoin(string(x0_21)),strjoin(string(x6)),16,m6,n6,'MS',f6}];


disp(['******************************************']);

% Test n°4
disp(['******************************************']);
disp(['Test n°4 : f2 et point réalisable']);

[x7,l7,m7,n7,f7] = LA(x0_22,lambda,mu,alpha,beta,eta_chap,tau,delta,@LA2,@gradLA2,@hessLA2,@c2,1);
res = [res;{'LA2',strjoin(string(x0_22)),strjoin(string(x7)),17,m7,n7,'Cauchy',f7}];

[x8,l8,m8,n8,f8] = LA(x0_22,lambda,mu,alpha,beta,eta_chap,tau,delta,@LA2,@gradLA2,@hessLA2,@c2,2);
res = [res;{'LA2',strjoin(string(x0_22)),strjoin(string(x8)),18,m8,n8,'MS',f8}];

disp(['******************************************']);

cres = cell2table(res);
cres.Properties.VariableNames = {'Fonction','Point_de_depart','Point_minimisant','lambda','mu','Nb_iterations','Pas_de_recherche','flag'};
writetable(cres,'testsLagrangienAugmente.xls','Sheet',1,'Range','A1');
disp(cres)