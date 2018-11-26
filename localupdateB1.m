function [B1] = localupdateB1(X, g1, hH, type)
B1 = [];

for j = g1 
    b1 = glmfit(hH, X(:,j),type{1},'link',type{2}, 'constant', 'off');
    B1 = [B1, b1];
end