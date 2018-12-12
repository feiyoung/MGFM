function [B1] = localupdateB1(X, g1, hH, type)

B1 = [];
n = size(X,1);
if strcmp(type{1}, 'binomial')    
for j = g1
    ntrail_j = length(unique(X(:,j)))-1;
    b1 = glmfit(hH, [X(:,j), ntrail_j*ones(n,1)],type{1},'link',type{2}, 'constant', 'off');
    B1 = [B1, b1];
end
else
    for j = g1
    b1 = glmfit(hH, X(:,j),type{1},'link',type{2}, 'constant', 'off');
    B1 = [B1, b1];
    end
end