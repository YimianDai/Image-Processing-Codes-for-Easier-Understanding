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

clc;
clear;
close all;

pr = 0.2; % percentage of r, pr = r / m
pe = 0.1; % percentage of sparse noise, #sparse noise = pe * m^2
m = 400;
r = round(pr * m);
A = randn(m, r);
B = randn(m, r);
X = A * B';
num_nz_E = round(pe * m^2); % number of non-zeros in E

E = zeros(m, m);
rand_ind = randperm(m^2);
E(rand_ind(1:num_nz_E)) = -5 + 10.*rand(num_nz_E, 1);

Y = X + E;
C = sqrt(m);
epsilon = eps();
lambda = 1 / sqrt(m);
[X_hat, E_hat] = wnnm_rpca_admm(Y, C, epsilon, lambda);

relative_error = sqrt(sum(sum((X_hat - X).^2))) / sqrt(sum(sum(X.^2)));
disp(['Relative Error: ' num2str(relative_error)]);

