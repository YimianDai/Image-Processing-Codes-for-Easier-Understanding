function X_hat = prox_reweighted_WNNP(Y, C, epsilon)
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

[U, S, V] = svd(Y, 'econ');
sigma_Y_arr= diag(S);
c1_arr = sigma_Y_arr - epsilon;
c2_arr = (sigma_Y_arr + epsilon).^2 - 4 * C;

idx = find(c2_arr > 0); 
psv_sigma_X_arr = (c1_arr(idx) + sqrt(c2_arr(idx))) / 2;

X_hat = U(:, idx) * diag(psv_sigma_X_arr) * V(:, idx)';

