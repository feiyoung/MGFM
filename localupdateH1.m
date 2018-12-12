function [hH] = localupdateH1(X, gcell, hB, type, dropout)
n = size(X, 1);
q = size(hB, 2);
ng = size(type,1);
if all(dropout ~=0) && (~isempty(setdiff(dropout , 1:ng)))
    error('dropout setting is wrong!')
end
% dlcell = cell(ng, 2);
% for j = 1:ng
%     [dlcell{j,1}, dlcell{j,2}] = typematchlink(type{j});
% end
Harray = zeros(n,q, ng);
w = ones(1,n);
for j = 1:ng
    H2 = [];
    if strcmp(type{j,1}, 'normal')
        w = 1./ (std(X(:, gcell{j})).^2);
        for i = 1:n
            h2 = glmfit(hB(gcell{j},:), X(i,gcell{j})',type{j,1},'link',type{j,2}, 'constant', 'off', 'weights', w);
            H2 = [H2, h2];
        end
    elseif strcmp(type{j,1}, 'binomial')
         p1 = length(gcell{j});
        for i = 1:n
            ntrail_i = length(unique(X(i,gcell{j})))-1;
           
            h2 = glmfit(hB(gcell{j},:), [X(i,gcell{j})',ntrail_i*ones(p1,1)],type{j,1},'link',type{j,2}, 'constant', 'off');
            H2 = [H2, h2];
        end
    else
        for i = 1:n
            h2 = glmfit(hB(gcell{j},:), X(i,gcell{j})',type{j,1},'link',type{j,2}, 'constant', 'off');
            H2 = [H2, h2];
        end
    end
    
    Harray(:,:,j) = H2';
end
idres = setdiff(1:ng, dropout);
hH = mean(Harray(:,:,idres), 3);