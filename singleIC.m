function q = singleIC(X_resamp, group, type, q_set, dropout, dc_eps, maxIter)

%dropout = 0; dc_eps = 1e-6;maxIter=50;
if(~exist('q_set', 'var') || isempty(q_set))
   q_set = 1:10;
end
p = size(X_resamp,2); omega = p^(-1/2);
q_num = length(q_set);
allhH = cell(q_num,1);
allhB = cell(q_num,1);
for r = 1:q_num
    %fprintf('r = %d \n ', r);
    [hH, hB] = gfm_evaluate(X_resamp, group, type, q_set(r), dropout, dc_eps, maxIter, omega);
    allhH{r} = hH;
    allhB{r} = hB;
end
Vr = zeros(2, q_num);
for r = 1:q_num
    %fprintf('r = %d \n ', r)
    Vr(:, r) = ICriteria(X_resamp, allhB{r}, allhH{r}, r, group, type);
end
IC = sum(Vr);
[~, id] = min(IC);
q = q_set(id);