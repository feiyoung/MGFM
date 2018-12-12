
%% ------------ Example 1. Generalized factor model: all are normal with homoskedasticity
clear;
i = 4; p = 100; n = 100; q = 6; % This is convergent! Curious.
[X, B, H] = gendat4(i,n,p);

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


%% ------------ Example 2. Generalized factor model: all are normal with  heteroskedasticity
clear;
i = 4; p = 100; n = 100; q = 6; % This is convergent! Curious.
[X, B, H] = gendat42(i,n,p);

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

[hH1, hB1] = factorm(X, 6);
[measurefun(H, hH), measurefun(H, hH1); ...
measurefun(B, hB), measurefun(B, hB1)]
%% ---------- Example 3. Generalized factor model: poisson + binary
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
%% ---------- Example 4. Generalized factor model: normal + poisson
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
%% ---------- Example 5. Generalized factor model: normal + poisson + binomial types
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
[measurefun(H, hH),measurefun(H, hH1); ...
measurefun(B, hB), measurefun(B, hB1)]

%% ---------- Example 6. Generalized factor model: pure binomial types
clear;
% generate data
p = 300; n = 300;
q = 6;
% generate H, B
rng(1); % since B is a fixed matrix.
ar_mat = cov_mat(ones(p,1)*sqrt(6), 0.5); % p*AR(1) covariance matrix
Z =  mvnrnd(zeros(1,p), ar_mat, n);
[Zdecomp,~] = eig(Z*Z');
B = sqrt(1/n)*Z'* Zdecomp(:,end:-1:end-q+1); % sort the eigenvectors by decreasing eigen values.
sB = sign(B(1,:));
B = B.* repmat(sB,p,1);  % ensure B is unique.
% generate H, B
rng(1);
H = mvnrnd(zeros(1,q),toeplitz(0.5.^(0:q-1)),n);
cF = cov(H, 0);
H = (H - repmat(mean(H),n,1))*cF^(-1/2);% ensure H is unqiue

% genarate X
g = 1:p;
mu = 1./(1 + exp(-H*B(g,:)')); % binary distribution.
N = 2;
X = binornd(N, mu);
%----gfm unified function
group = [ones(1,p)]; % full-one vector indicates all variables belong to the same type.
type = cell(1,2);
type{1,1} = 'binomial';  
type{1,2} = 'logit'; 
q = 6; % q is given
[hH, hB, history] = gfm(X, group, type, 6, 0, 1e-3, 6);
[hH1, hB1] = factorm(X, 6);
[measurefun(H, hH),measurefun(H, hH1); ...
measurefun(B, hB), measurefun(B, hB1)]
% q is estimated
[hH, hB, history] = gfm(X, group, type, [], 0, 1e-3, 6);
[measurefun(H, hH),measurefun(H, hH1); ...
measurefun(B, hB), measurefun(B, hB1)]