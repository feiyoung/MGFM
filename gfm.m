function [hH, hB, history] = gfm(X, group, type, q, dropout, dc_eps, maxIter, omega,N, q_set, output)
%% This function is used to conduct the Generalized Factor Model.
% Author: Liu Wei(117020208005@2017.swufe.edu.cn)
% Version: 2018-11-25
% Iterative algorithm to estimate the factors H and loading matrix B in
% generalized factor model.
% --X: a matrix with dimension of n*p(p=(p1+p2+..+p_d)), observational mixed data
% matrix, d is the types of variables, p_j is the dimension of j-th type
% variable.
% --group: a vector with length equal to p, specify each column of X belonging to which group.
% --type: a d-by-2 cell variable, specify the type of variable and link function on each
% group. For example, type={'poisson', 'log'; 'binomial', 'probit'}, and it
% is referred to the help file of glmfit() function for more details.
% --q: a positive integer or empty, specify the number of factors. If q is
% empty, then IC critera is used to dertemined q automatically.
% --dropout: a proper subset of [1, 2, ..., d],  specify which group to be
% dropped in obtaining the initial estimate of factor matrix H, and the aim
% is to ensure the convergence of algorithm leaded by weak signal variable
% types. Optional parameter with default as 0, no group dropping.
% -- maxIter: a positive integer, specify the times of iteration. Optional
% parameter with default as 50.
% -- dc_eps: a positive real number, specify the tolerance of varing quantity of objective
% function in the algorithm. Optional parameter with default as 1e-6.
% -- omega: a positive integer, correlated with p, to control the order of objective function.
% It does not have essential effect on the algorithm. Optional parameter
% with default as p^{-1/2}.
% -- N: a positive integer, specify the times of leave-one-out bootstrap if
% q is empty. Optional parameter with default as 1. 
% -- q_set: a positive integer set, specify the candidates of factor number
% q, (optional) default as [1:8] according Bai,2013.
% -- output: a logical value with true or false, specify whether ouput the
% mediate information in iteration process, (optional) default as true.
if(~exist('dropout', 'var') || isempty(dropout))
    dropout=0;
end
if(~exist('dc_eps', 'var')|| isempty(dc_eps))
    dc_eps = 1e-6; 
end
if(~exist('maxIter', 'var')|| isempty(maxIter))
    maxIter = 50; 
end
[n, p] = size(X);
if(~exist('omega', 'var')|| isempty(omega))
    omega = p^(-1/2);
end
if(~exist('N', 'var') || isempty(N))
    N=1;
end
if(~exist('q', 'var') || isempty(q))
    if(~exist('q_set', 'var') || isempty(q_set))
       q_set = 5:7;
    end
    q_vec = zeros(1, N);
    rng(1);
    fprintf('Start to determine the factor number q ....\n')
    for k = 1:N
        X_resamp = datasample(X, n-1, 'Replace',false);
        q_vec(k) = singleIC(X_resamp, group, type, q_set, dropout, dc_eps, maxIter);
    end
    q= round(mean(q_vec));
    fprintf('The factor number q is estimated as %d . \n', q);
end



if(~exist('output', 'var') || isempty(output))
    output = false;
end
fprintf('Start booting the iterative algorithm ....\n')
[hH, hB, history] = gfm_evaluate(X, group, type, q, dropout, dc_eps, maxIter, omega, output);
history.dc_eps = dc_eps; history.dropout = dropout; history.q = q;
end
