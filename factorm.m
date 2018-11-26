function [hatH, hatB] = factorm(X, q)
n = size(X, 1);
[evec,~] = eig(X*X');
hatH = evec(:, end:-1:end-q+1)* sqrt(n);
hatB = X'*hatH/n; % ensure B'B is a diagonal matrix with decreasing order 
sB = sign(hatB(1,:));
hatH = hatH.* repmat(sB,n,1); % ensure F is unqiue (A2) condition
cF = cov(hatH, 0);
hatH = (hatH - repmat(mean(hatH),n,1))*cF^(-1/2);
end

