%%% TV demo %%
%%% By Guy Gilboa
% http://visl.technion.ac.il/~gilboa/PDE-filt/tv_denoising.html

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on: [ROF92] L. Rudin, S. Osher, E. Fatemi,
% "Nonlinear Total Variation based noise removal algorithms",
% Physica D 60 259-268,1992.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Only add some comments for easier understanding
% It should be noted that the \lambda and u of Prof. Gilboa's code here is not  
% the same as the equations (2.9c), (2.6), (2.5a) and (2.8a)in the paper. 
% I wrote a blog (in Chinese) about it to express the reason, please refer to
% http://dym.mobi/post/research/rof92-tv-3

clc;
clear all;
close all;

addpath('support');
addpath('../../images');
addpath('../../utilities');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Simple TV denoising  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I = double(imread('cameraman.tif')); % load image
% I=double(I(20:120,10:105)); % cut a piece, convert to double

% params
iter = 80; 
% dt=0.2; eps=1; 

%%% Add noise
std_n = 20; % Gaussian noise standard deviation
In = randn(size(I))*std_n; % White Gaussian noise
I0 = I + In;  % noisy input image
% show original and noisy images
figure(1); imshow(uint8(I)); title('Original')
figure(2); imshow(uint8(I0)); title('Noisy image')
% denoise image by using tv for some iterations
J = tv(I0,iter); 
psnr1 = imgpsnr(I, J);
figure(3); imshow(uint8(J)); 
title(['Denoised image, PSNR = ' num2str(psnr1)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  TV denoising with (scalar) data fidelity term  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Here we assume the noise vaiance is known.
%% Therefore we can compute lambda (instead of fixing
%% iter). TV therfore can run automatically quite well.
%% The process is slower.
J = I0; 
% params
ep_J = 0.01; % minimum mean change in image J
lam = 0; % \lambda
J_old = 0;
iter = 10; 
dt = 0.2; 
eps = 1; % avoid nan (dividing zero)
var_n = std_n^2; % noise variance
lam_arr = [];
while (mean(mean(abs(J - J_old))) > ep_J),  % iterate until convergence
   J_old = J;
   J = tv(J,iter,dt,eps,lam,I0);     % scalar lam
   lam = calc_lam(J,I0,var_n,eps); % update lambda (fidelity term)
   lam_arr = [lam_arr lam];
end
psnr2 = imgpsnr(I, J);
figure(4); imshow(uint8(J));
title(['Denoised image with lambda' num2str(psnr2)]);
figure; plot(lam_arr);