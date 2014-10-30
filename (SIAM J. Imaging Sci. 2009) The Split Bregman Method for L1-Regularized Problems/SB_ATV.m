function u = SB_ATV(f,mu)
% Split Bregman Anisotropic Total Variation Denoising
%
% u = arg min_u |\nabla_x u| + |\nabla_y u| + mu/2||u-f||_2^2
%
% Refs:
%  *Goldstein and Osher, The split Bregman method for L1 regularized problems
%   SIAM Journal on Imaging Sciences 2(2) 2009
% 
% Previous version (original branch) of this algorithm code 
% Split Bregman method for Total Variation Denoising 
% by Benjamin Tremoulheac 
% http://www.mathworks.com/matlabcentral/fileexchange/authors/94260
%
% Motivation to fork Benjamin's codes:
% 1. there are several places that Benjamin's codes are not consistent with the equations
%    in the paper, for example the most important equation (4.2), which makes me feel 
%    confused when reads his codes.
% 2. Benjamin's codes are not easy to understand for lacking comments.
%
% Changes made by this version: 
% 1. Add detailed comments for easier understanding
% 2. Change the variable names to be consistent with the symbols in the paper
% 3. Rewrite some codes to be consistent with the algorithm steps in the paper
%
% More information about this version (in Chinese), please refer to 
% 
% 
% INPUT ARGUMENTS : f  - the noisy image (gray-level scale)
%                   mu - .
%
% OUTPUT ARGUMENTS : u - denoised image.

f = f(:); % reshape the noisy image into a column vector 
          % the reason why we need a column vector is that 
          % the step u^{k+1} = G_{i,j}^{k} in the paper is solved by Gauss-Seidel solution 
          % but here we use the cgs function - Conjugate gradients squared method by MATLAB
          % cgs solve the system of linear equations A*X=B, X is u^{k+1}
          % because in cgs function, B must be a column vector, so u is also a column vector
          % similarily, b_x^{k}, b_y^{k}, d_x^{k}, d_y^{k, s^{k} are also column vectors 
n = length(f); % n is the length of one side of the image
               % here, we implicitly assume the image is squared (n*n pixels)
[nablaX nablaY laplace] = diffoper(sqrt(n));
% initialize
bX = zeros(n,1);
bY = zeros(n,1);
dX = zeros(n,1);
dY = zeros(n,1);
u  = f;
err = 1;
tol = 1e-3;
lambda = 0.05;
k = 1; 
while err > tol
    fprintf('it. %g ',k);
    up = u; % up is the u of previous iteration
    
    % update u^{k+1}
    % solve the system of linear equations A*X=B 
    cgsA = mu*speye(n) - lambda*laplace;
    cgsB = mu*f + lambda*nablaX'*(dX-bX) + lambda*nablaY'*(dY-bY);
    [u,~] = cgs(cgsA, cgsB, 1e-3, 100); % X = cgs(A,B,TOL,MAXIT)      
    
    dX = shrink(nablaX*u + bX, 1/lambda);
    
    dY = shrink(nablaY*u + bY, 1/lambda);
    
    % update b_x^{k+1}
    bX = bX + (nablaX*u - dX);
    
    % update b_y^{k+1}
    bY = bY + (nablaY*u - dY);    

    err = norm(up-u)/norm(u);
    fprintf('err=%g \n',err);
    k = k+1;
end
fprintf('Stopped because norm(up-u)/norm(u) <= tol=%.1e\n',tol);
end

function [nablaX nablaY laplace] = diffoper(N)
% INPUT ARGUMENTS : N - the side length of a squared image
%
% OUTPUT ARGUMENTS : 
% nablaX - the \nabla_x operator matrix to get the difference result in the X direction
% nablaY - the \nabla_y operator matrix to get the difference result in the Y direction
% laplace - the Laplace operator matrix (\Delta)
% 
% Note that: 
% those difference operators are defined in the matrix form,
% not the explict rule we often used : latterItem - formerItem
%
% Tips: 
% since difference operators are also linear operators 
% and linear operators can be written as a matrix
D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
D(:,1) = []; D(1,1) = 0;
nablaX = kron(speye(N),D);
nablaY = kron(D,speye(N));
laplace = -[nablaX; nablaY]' * [nablaX; nablaY]; % it can also be written as follows
% laplace = nablaX'*nablaX + nablaY'*nablaY; 
end  

function output = shrink(x,gamma)
% equation (3.10) in the paper 
output = max(abs(x)-gamma, 0) .* sign(x);
end