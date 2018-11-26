
%% ------------ Example 1. Generalized factor model: all are normal
clear;
i = 4; p = 100; n = 100; q = 6; 
[X, B, H] = gendat4(1,i,n,p);

% unified function
group = [ones(1,p)]; % full-one vector indicates all variables belong to the same type.
type = cell(1,2);
type{1,1} = 'normal';  % the type is 'normal'
type{1,2} = 'identity'; % the link funciton is 'identity'
q= 6; 
[hH, hB, history] = gfm(X, group, type, q); % start to estimate.
%
q = []; % use IC criteria with leave-one-out bootstrap to dertermine q
[hH, hB, history] = gfm(X, group, type, q); % start to estimate.
% compare the precision of estiamtion
[measurefun(H, hH),measurefun(H, hH1); ... 
measurefun(B, hB), measurefun(B, hB1)]
%% ---------- Example 2. Generalized factor model: poisson + binary
clear;
i = 1; p = 100; n = 100;
[X, H, B,hB] = gendat(i, n, p);
measurefun(B, hB)
[hH1, hB1] = factorm(X, 6);
hH1(1:6, 1:6)
H(1:6, 1:6)

% unified function
type = cell(2,2);
type{1,1} = 'poisson'; type{1,2} = 'log';
type{2,1} = 'binomial';  type{2,2} = 'probit'; % Also can choose 'logit'
group = [ones(1,floor(p/2)), 2*ones(1, p-floor(p/2))]; % The first p/2 variables are of 'poisson'
% type, and the last p/2 variables are of 'binomial' type.
dropout=[2]; % Because 'binomial' type is weak signal relative to 'poisson' type, we add the 'dropout'
% step!
maxIter = 10;
q= 6; 
[hH, hB, history] = gfm(X, group, type, q, dropout, [], maxIter);
history
[measurefun(H, hH),measurefun(H, hH1); ...
measurefun(B, hB), measurefun(B, hB1)]

q = []; % use IC criteria with leave-one-out bootstrap to dertermine q
[hH, hB, history] = gfm(X, group, type, q, dropout); %%
% compare the precision of estiamtion
[measurefun(H, hH),measurefun(H, hH1); ... 
measurefun(B, hB), measurefun(B, hB1)]
%% ---------- Example 3. Generalized factor model: normal + poisson
clear;
i =3; p = 50; n = 50;
[X, H, B,hB] = gendat6(i, n, p);
[hH1, hB1] = factorm(X, 6);
measurefun(H, hH1)
measurefun(B, hB1)
% unified function, well done, it is very very good.
group = [ones(1,floor(p/2)), 2*ones(1, p-floor(p/2))];
type = {'normal', 'identity'; 'poisson', 'log'}; q= 6;
[hH, hB, history] = gfm(X, group, type, q);
history
[measurefun(H, hH),measurefun(H, hH1); ...
measurefun(B, hB), measurefun(B, hB1)]

q = []; % use IC criteria with leave-one-out bootstrap to dertermine q
[hH, hB, history] = gfm(X, group, type, q); %%
% compare the precision of estiamtion
[measurefun(H, hH),measurefun(H, hH1); ... 
measurefun(B, hB), measurefun(B, hB1)]
%% ---------- Example 4. Generalized factor model: normal + poisson + binomial types
clear;
n = 200; p=200; i = 1;
[X, H, B,hB] = gendat7(i, n, p);
[hH1, hB1] = factorm(X, 6);
measurefun(H, hH1)
measurefun(B, hB1)
% unified functions test
type = cell(3,2);
type{1,1} = 'normal'; type{1,2} = 'identity';
type{2,1} = 'poisson'; type{2,2} = 'log';
type{3,1} = 'binomial';  type{3,2} = 'logit';
group = [ones(1,floor(p/3)), 2*ones(1, floor(2*p/3)-floor(p/3)), 3*ones(1, p-floor(2*p/3))];
q= 6; dropout=[3];omega = p^(-1);
[hH, hB, history] = gfm(X, group, type, q, dropout, [], [], omega);
history
[measurefun(H, hH),measurefun(H, hH1); ...
measurefun(B, hB), measurefun(B, hB1)]

q = []; % use IC criteria with leave-one-out bootstrap to dertermine q
[hH, hB, history] = gfm(X, group, type, q, dropout); %%
% compare the precision of estiamtion
[measurefun(H, hH),measurefun(H, hH1); ... 
measurefun(B, hB), measurefun(B, hB1)]