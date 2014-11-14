% This file is an implementation of ROF92 TV denoising algorithm. 
% 
% Refs:
%  *Rudin L I, Osher S, Fatemi E. 
%   Nonlinear total variation based noise removal algorithms[J]. 
%   Physica D: Nonlinear Phenomena, 1992, 60(1): 259-268.
%
% Note that according to whether the derivation of equation (2.5a) in the 
% continuous form is done completely, I write two versions of ROF denoising
% code. This version (version 1) directly discretize equation (2.5a).
%
% For more information and details about this algorithm, please refer to 
% http://dym.mobi/post/research/rof92-tv-v1
%
% I'm bound to say that the denoising performance of this code version is not
% very good and I can't find the reason now. I still upload this code in the
% hope that someone would tell me where is the problem lying in. 
% 
% If you have some idea to improve this code, please let me know.
% My email is daiyimian@outlook.com
%
% For better denosing performance, you can refer to my other TV codes such as 
% (SIAM J. Imaging Sci. 2009) 
% The Split Bregman Method for L1-Regularized Problems
% or the version 2.

clc;
clear all;
close all;

addpath('support');
addpath('../../images');
addpath('../../utilities');

img = double(imread('cameraman','tif')); % noise-free image
sigma = 20; % standard variance of the AWGN 
nimg = img + sigma*randn(size(img)); % noisy image
u0 = nimg; % noisy image, the same notation as it in the paper
u = u0; % iterated image u in the paper
dt = 0.2; % \Delta t in the paper, time step
h = 1; % h in the paper
epsilon = 10^-3; % \epsilon to avoid nan (not a number) occurrence
                 % nan occurs when dividing zero
iterNum = 100;                 
lambda_arr = [];
for i = 1:iterNum
    % update lambda according to equation (2.9c) in the paper 
    lambda = complambda(u,u0,h,sigma);  
    lambda_arr = [lambda_arr lambda]; % record every lambda to see its changes
    
    % the denominator of first large fraction in (2.8a)    
    den1  = sqrt(backdiffx(u).^2 + minmod(backdiffy(u),frontdiffy(u)).^2 + ...
                 epsilon^2); 
    % the denominator of second large fraction in (2.8a)
    den2  = sqrt(backdiffy(u).^2 + minmod(backdiffx(u), frontdiffx(u)).^2 + ...
                 epsilon^2);
    frac1 = backdiffx(u) ./ den1; % the first large fraction in (2.8a)
    frac2 = backdiffy(u) ./ den2; % the second large fraction in (2.8a)
    % update u according to equation (2.8a) in the paper
    uNext = u + dt/h*(frontdiffx(frac1) + frontdiffy(frac2)) - dt*lambda*(u-u0);
    u = uNext;
end
denoImg = u;

nimgPSNR = imgpsnr(img,nimg);
denoPSNR = imgpsnr(img,denoImg);

subplot(1,3,1); imshow(uint8(img)); 
title('Noise-free image');
subplot(1,3,2); imshow(uint8(nimg)); 
title(['Noisy image, PSNR = ' num2str(nimgPSNR)]);
subplot(1,3,3); imshow(uint8(denoImg)); 
title(['denoised image, PSNR = ' num2str(denoPSNR)]);

figure; plot(lambda_arr);