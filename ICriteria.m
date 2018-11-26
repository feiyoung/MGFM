function Vr = ICriteria(X, hB, hH, r, group, type)

[n, p] = size(X);
omega = 1/p;
ind_set = unique(group);
ng = length(ind_set);
gcell = cell(1, ng);
for j = 1:ng
    gcell{j} = find(group==j);
end
c = objfunc(hH, hB, X, omega, gcell, type);
Vr = [log(c) , r/min(sqrt(n), sqrt(p))^2*log(min(sqrt(n), sqrt(p))^2)];
%Vr = [log(c) , r*(n+p)/(n*p)*log(min(sqrt(n), sqrt(p))^2)];
%Vr = [log(c) , r*(n+p)/(n*p)*log(n*p/(n+p))];
%Vr = [c , r*0.6*(n+p)/(n*p)*log(n*p/(n+p))];