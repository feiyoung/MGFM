function [X, H, B,hB] = gendat7(i, n, p)
% This function generate simulated data with sample size n and dimension p.
% three exponential family distribution: Normal + Poisson + binary
if(~exist('p', 'var'))
   p = 50; 
end
if(~exist('n', 'var'))
   n = 300; 
end
q = 6;
d = 2;
% generate H, B
rng(1); % since B is a fixed matrix.
ar_mat = cov_mat(ones(p,1)*sqrt(4), 0.5); % p*AR(1) covariance matrix
Z =  mvnrnd(zeros(1,p), ar_mat, n);
[Zdecomp,~] = eig(Z*Z');
B = sqrt(1/n)*Z'* Zdecomp(:,end:-1:end-q+1); % sort the eigenvectors by decreasing eigen values.
sB = sign(B(1,:));
B = B.* repmat(sB,p,1);  % ensure B is unique.
% generate H, B
rng(i);
H = mvnrnd(zeros(1,q),toeplitz(0.5.^(0:q-1)),n);
cF = cov(H, 0);
H = (H - repmat(mean(H),n,1))*cF^(-1/2);% ensure H is unqiue

% genarate X
g1 = 1:floor(p/3); % identity, normal
g2 = (floor(p/3)+1):floor(2*p/3); % Poisson: exp
g3 = (floor(2*p/3) + 1):p; % Bernoulli

mu1 = H* B(g1,:)'; % normal dstribution.
X1 = normrnd(mu1, 1); 
mu2 = exp(H* B(g2,:)'); % poisson distribution
X2 = poissrnd(mu2);
mu3 = 1./(1 + exp(-H*B(g3,:)')); % binary distribution.
X3 = binornd(1, mu3);

X = [X1 X2 X3];
B1 = [];
for j = g1
b1 = glmfit(H, X(:,j),'normal', 'constant', 'off');
B1 = [B1, b1];
end

B2 = [];
for j =g2
b2 = glmfit(H, X(:,j),'poisson','link','log', 'constant', 'off');

B2 = [B2, b2];
end

B3 = [];
for j =g3

b3 = glmfit(H, X(:,j), 'binomial', 'link', 'logit', 'constant', 'off');
B3 = [B3, b3];
end

hB = [B1,B2, B3]';
end