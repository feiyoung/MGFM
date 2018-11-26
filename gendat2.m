function [X, B, H] = gendat2(Case,seed,n,p)
%% linear factor model to generate data.

q = 6;
d = 2;
if(~exist('n', 'var'))
    n = 300;
end
if(~exist('p', 'var'))
    p = 50;
end
%-------------------- Case I -----------------------%
rng(seed);  % For reproducibility
ar_mat = p*toeplitz(0.5.^(0:p-1)); % p*AR(1) covariance matrix
Z =  mvnrnd(zeros(1,p), ar_mat, n);
[Zdecomp,~] = eig(Z*Z');
% Zdecomp' * Zdecomp
B = sqrt(1/n)*Z'* Zdecomp(:,end:-1:end-q+1); % sort the eigenvectors by decreasing eigen values.
switch Case
    case 1
        H = mvnrnd(zeros(1,q),toeplitz(0.5.^(0:q-1)),n);
        sB = sign(B(1,:));
        B = B.* repmat(sB,p,1); 
        cF = cov(H, 0);
        H = (H - repmat(mean(H),n,1))*cF^(-1/2);% ensure F is unqiue (A2) condition
        
        X = H * B' + mvnrnd(zeros(1,p), diag(ones(p,1)), n);
       
    case 2
        %----------------------------Case II --------------------------%
        H12 =  mvnrnd(zeros(1,2), toeplitz(0.5.^(0:1)), n);
        H3 = abs(sum(H12,2)) + H12(:,1) .* randn(n,1);
        H4 = (sum(H12,2)).^2 + H12(:,2) .* randn(n,1);
        H5 = binornd(1, exp(H12(:,2))./ (1 + exp(H12(:,2))));
        H6 = binornd(1, normcdf(H12(:,2)));
        H = [H12 H3 H4 H5 H6];
        sB = sign(B(1,:));
        B = B.* repmat(sB,p,1); 
        cF = cov(H, 0);
        H = (H - repmat(mean(H),n,1))*cF^(-1/2);% ensure F is unqiue (A2) condition
        
        X = H * B' + mvnrnd(zeros(1,p), diag(ones(p,1)), n);
        
end
end