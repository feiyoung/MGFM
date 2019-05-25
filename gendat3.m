function [Y, X, H, B,hB] = gendat3(i, p, n)
% This function generate simulated data with sample size n and dimension p.
% two exponential family distribution: Normal and binary.
%----------------------------------------------------------------------------------
% Author:       Liu Wei
% Maintainer:    Liu Wei <weiliu@smail.swufe.edu.cn>
% Date:         Jan. 13, 2018
% Copyright (c) 2018, Liu Wei
% All rights reserved.
if(~exist('p', 'var'))
   p = 50; 
end
if(~exist('n', 'var'))
   n = 300; 
end
q = 6;
d = 2;
% generate H, B
rng(i);
ar_mat = cov_mat(ones(p,1)*sqrt(p), 0.5); % p*AR(1) covariance matrix
Z =  mvnrnd(zeros(1,p), ar_mat, n);
[Zdecomp,~] = eig(Z*Z');
B = sqrt(1/n)*Z'* Zdecomp(:,end:-1:end-q+1); % sort the eigenvectors by decreasing eigen values.
H = mvnrnd(zeros(1,q),toeplitz(0.5.^(0:q-1)),n);
sB = sign(B(1,:));
B = B.* repmat(sB,p,1);  % ensure B is unique.
cF = cov(H, 0);
H = (H - repmat(mean(H),n,1))*cF^(-1/2);% ensure H is unqiue

% genarate X
g2 = 1:p; 
mu2 = 1./(1 + exp(-H*B(g2,:)')); % binary distribution.
X2 = binornd(1, mu2);
X = X2;

B2 = [];
for j =g2
%b2 = glmfit(H, X(:,j),'poisson','link','log', 'constant', 'off');
b2 = glmfit(H, X(:,j), 'binomial', 'link', 'logit', 'constant', 'off');
B2 = [B2, b2];
end
hB = B2';

% generate beta, W, epsilon, H(.), G(.), Y
beta = [1,0,1,1,1,1;0,1,-1,1,-1,1];
W = H* beta';
epsi = randn(n,1);
Y = 2*W(:,1).^2 + 4*abs(W(:,2) + 1) + W(:,2).^2 .* epsi;  % 

end