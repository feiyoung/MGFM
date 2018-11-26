function [X, H, B,hB] = gendat5(i, n, p)
% This function generate simulated data with sample size n and dimension p.
% two exponential family distribution: possion
if(~exist('p', 'var'))
   p = 50; 
end
if(~exist('n', 'var'))
   n = 300; 
end
q = 6;
d = 2;
rng(1); % since B is a fixed matrix.
ar_mat = cov_mat(ones(p,1)*sqrt(1), 0.5); % p*AR(1) covariance matrix
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
g2 = 1:p; 
mu3 = exp(H* B(g2,:)'); % poisson distribution
X3 = poissrnd(mu3);
B3 = [];
for j =g2
b3 = glmfit(H, X3(:,j),'poisson','link','log', 'constant', 'off');
B3 = [B3, b3];
end
hB = B3';
X = X3;
end