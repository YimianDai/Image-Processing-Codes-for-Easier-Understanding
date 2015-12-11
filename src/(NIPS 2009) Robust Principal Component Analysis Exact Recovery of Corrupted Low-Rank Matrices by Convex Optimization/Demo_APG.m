%% Create artificial low-rank image
perHei = 10;
perWid = 20;
hei = 50*perHei;
wid1 = 5*perWid;  rect1 = ones(hei, wid1) * 75;
wid2 = 15*perWid; rect2 = ones(hei, wid2) * 184;
wid3 = 8*perWid;  rect3 = ones(hei, wid3) * 223;
wid4 = 18*perWid; rect4 = ones(hei, wid4) * 113;
wid5 = 10*perWid; rect5 = ones(hei, wid5) * 38;
A = uint8([rect1 rect2 rect3 rect4 rect5]);
D = imnoise(A, 'salt & pepper');
A = double(A);
D = double(D);
E = D - A;
%% Create low-rank matrix
% m = 100; % m = 100, 200, 400, 800
% r = 5;  % rank(A) = 5, 10, 20, 40
% U = randn(m,r);
% V = randn(m,r);
% A = U*V';
% E = (imnoise(zeros(size(A)),'salt & pepper',0.1) > 0) .* (rand(size(A))-0.5)*2*500;
% D = A + E;
%% Apply APG to recover A_hat and E_hat
[rows, cols] = size(D);
lambda = rows^(-1/2);
tic
[A_hat,E_hat,numIter] = proximal_gradient_rpca(D, lambda);
toc

norm(A_hat-A, 'fro') / norm(A, 'fro')

figure; imshow(A,[]); title('original A');
figure; imshow(A_hat,[]); title('recovered A');