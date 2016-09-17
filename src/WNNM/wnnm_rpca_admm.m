function [X, E] = wnnm_rpca_admm(Y, C, epsilon, lambda, tol, max_iter)
% \min_{X, E} \lambda \Vert E \Vert_1 + \Vert X \Vert_{w,*}, s.t. Y = X + E
% Reference:
% [1] Gu, S.; Xie, Q.; Meng, D.; Zuo, W.; Feng, X. & Zhang, L. 
%     Weighted Nuclear Norm Minimization and Its Applications to Low Level Vision 
%     International Journal of Computer Vision, 2016, 1-26
% If you find any bug or error, please pull a issue or email me 
% (yimian.dai at gmail dot com).
% Project page: 
% https://github.com/YimianDai/Image-Processing-Codes-for-Easier-Understanding
% This code is my personal practice code. The code of the authors of Ref. [1] 
% could be found in http://www4.comp.polyu.edu.hk/~cslzhang/papers.htm.

% Initialize 
[m, n] = size(Y);
if ~exist('max_iter', 'var')
    max_iter = 1000;
end
if ~exist('tol', 'var')
    tol = 10.0^(-7);
end
if ~exist('lambda', 'var')
    lambda = 1 / sqrt(m);
end
if ~exist('epsilon', 'var')
    epsilon = eps();
end
if ~exist('C', 'var')
    C = sqrt(m);
end

norm_fro = norm(Y, 'fro');
norm_two = norm(Y, 2); 

X = Y;
rho = 1.05;
mu = 1 / norm_two;
L = zeros(m, n);
% norm_inf = norm(Y(:), Inf);
% dual_norm = max(norm_two, norm_inf);
% L = Y / dual_norm;
iter = 0;
converged = false;

while ~converged
    iter = iter + 1;    
    
    E = prox_l1(Y + L / mu - X, lambda / mu);    
    
    X = prox_reweighted_WNNP(Y + L / mu - E, C / mu, epsilon);    
    
    L = L + mu * (Y - X - E);    
    
    mu = rho * mu;
    
    stopCriterion = norm(Y - X - E, 'fro') / norm_fro;
    if stopCriterion < tol
        converged = true;
    end
    
    if mod(iter, 10) == 0
        disp(['#Iteration ' num2str(iter) ' r(X) ' num2str(rank(X)) ...
            ' |E|_0 ' num2str(length(find(abs(E)>0))) ...
            ' stopCriterion ' num2str(stopCriterion)]);
    end  
    
    if ~converged && iter >= max_iter
        disp('Maximum iterations reached') ;
        converged = true;       
    end    
end
    


    