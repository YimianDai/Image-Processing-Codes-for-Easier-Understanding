% This file demonstrates the Split Bregman method for Total Variation denoising
%
% SB_ATV.m  Split Bregman Anisotropic Total Variation Denoising
% SB_ITV.m  Split Bregman Isotropic Total Variation Denoising
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

clc; 
clear all;
close all;

addpath('../../images');
addpath('../../matlab_support');

img = double(imread('Lena512','png'));
sigma = 20;
nimg = img + sigma*randn(size(img));

mu = 0.05;

atvImg = SB_ATV(nimg,mu);
atvImg = reshape(atvImg, size(img));
itvImg = SB_ITV(nimg,mu);
itvImg = reshape(itvImg, size(img));

% fprintf('ATV Rel.Err = %g\n',norm(g_denoise_atv(:) - img(:)) / norm(img(:)));
% fprintf('ITV Rel.Err = %g\n',norm(isoTVImg(:) - img(:)) / norm(img(:)));

figure; colormap gray;
subplot(221); imagesc(img); axis image; title('Original');
subplot(222); imagesc(nimg); axis image; title('Noisy');
subplot(223); imagesc(atvImg); axis image; 
title('Anisotropic TV denoising');
subplot(224); imagesc(itvImg); axis image; 
title('Isotropic TV denoising');

imgpsnr(img, atvImg)
imgpsnr(img, itvImg)