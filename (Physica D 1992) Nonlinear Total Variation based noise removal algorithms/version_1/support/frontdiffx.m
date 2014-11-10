function frontDiffUx = frontdiffx(u)
% forward diff in x direction in equation (2.9a)
% \Delta^x_{-}u_{ij} = -(u_{i-1,j} - u_{ij}) 

rows = size(u,1);
% u([1 1:rows-1], :) is u_{i-1,j} in (2.9a)
frontDiffUx = -(u([1 1:rows-1], :) - u); 
