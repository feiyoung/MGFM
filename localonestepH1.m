function H2 = localonestepH1(X, B, H, gcell, type)
% update H
% nonlinear factor model: poisson + binary
[n, p] = size(X);
q = size(H, 2);
ng = size(type,1);
mucell = cell(1,ng);
for j = 1:ng
    switch type{j,1}
    case 'normal'
        mucell{j} = H * B(gcell{j},:)';
    case 'poisson'
        mucell{j} = exp(H * B(gcell{j},:)');
    case 'binomial'
        Xj = X(:, gcell{j});
        ntrail_j = length(unique(Xj(:,1)))-1;
        mucell{j} =ntrail_j*1./(1 + exp(-H * B(gcell{j},:)'));
    end
end
% % the result is right
d2f = cell(n,1);
for i = 1:n
    Bng = zeros(q,q);
    for j = 1:ng
        switch type{j,1}
    case 'normal'
        %W = diag(1./ (std(X(:,gcell{j})).^2));
        Bng = Bng + B(gcell{j},:)'*B(gcell{j},:);
    case 'poisson'
        Bng = Bng + (repmat(mucell{j}(i,:), q, 1)'.* B(gcell{j},:))'* B(gcell{j},:);
    case 'binomial'
        Bng = Bng + (repmat(mucell{j}(i,:), q,1)' .* B(gcell{j},:))' *(repmat(1-mucell{j}(i,:), q,1)' .* B(gcell{j},:));
        end
    end
    d2f{i} = Bng;
end
df2 = zeros(n, q);
for j = 1:ng
    df2 = df2 + (X(:, gcell{j})- mucell{j})* B(gcell{j},:);
end
cell1 = cell(n,1);
for i = 1:n
    cell1{i} = df2(i,:)';
end
H2update = cellfun(@(x, y) x\y, d2f, cell1, 'UniformOutput', 0);
% H2 = H - cell2array(H2update)';
H2updatemat = reshape(cell2mat(H2update), q,n);
H2 = H - H2updatemat';
