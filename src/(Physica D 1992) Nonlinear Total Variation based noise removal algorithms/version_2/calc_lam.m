function lam=calc_lam(I,I0,var_n,ep)
% Only add some comments to Prof. Gilboa's code.
%
% private function: calc_lam (by Guy Gilboa).
% calculate scalar fidelity term of the Total-Variation
% input: current image, original noisy image, noise variance, epsilon
% output: lambda
% example: lam=calc_lam(I,I0,var_n,ep)

if ~exist('ep') % good for 256 gray-level
   ep=1;
end

[ny,nx] = size(I); 
ep2 = ep^2;

%% estimate derivatives
% central differencing, I_x = \frac{I_{i,j+1} - I_{i,j-1}}{2}
I_x = (I(:,[2:nx nx])-I(:,[1 1:nx-1]))/2; 
% central differencing, I_y = \frac{I_{i+1,j} - I_{i-1,j}}{2}
I_y = (I([2:ny ny],:)-I([1 1:ny-1],:))/2;
% I_{xx} = \frac{I_{i,j+1}+I_{i,j-1}-2I}{2}
I_xx = I(:,[2:nx nx])+I(:,[1 1:nx-1])-2*I;
% I_{yy} = \frac{I_{i+1,j}+I_{i+1,j}-2I}{2}
I_yy = I([2:ny ny],:)+I([1 1:ny-1],:)-2*I;
% I_{xy} = \frac{I_{i+1,j+1}+I_{i-1,j-1}-I_{i+1,j}-I_{i,j+1}}{4}
Dp = I([2:ny ny],[2:nx nx])+I([1 1:ny-1],[1 1:nx-1]);
Dm = I([1 1:ny-1],[2:nx nx])+I([2:ny ny],[1 1:nx-1]);
I_xy = (Dp-Dm)/4;

%% compute numerator and denomenator
Num = I_xx.*(ep2+I_y.^2)-2*I_x.*I_y.*I_xy+I_yy.*(ep2+I_x.^2);
Den = (ep2+I_x.^2+I_y.^2).^(3/2);
Div = Num./Den;

%% fidelity term
% It should be noted that the \lambda of Prof. Gilboa's code here is not  
% the same as the equation (2.6) and (2.9c). 
% (2.9c) is the discrete form of (2.6).
% I wrote a blog (in Chinese) about it to express the reason, please refer to
% http://dym.mobi/post/research/rof92-tv-3

lam = mean(mean(Div.*(I-I0)))./var_n; 
 