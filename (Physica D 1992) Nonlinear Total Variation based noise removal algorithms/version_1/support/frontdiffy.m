function frontDiffUy = frontdiffy(u)
% forward diff in y direction similar to equation (2.9a)
% \Delta^y_{-}u_{ij} = -(u_{i,j-1} - u_{ij}) 

cols = size(u,2);
% u(:, [1 1:cols-1]) is u_{i,j-1} as (2.9a)
frontDiffUy = -(u(:, [1 1:cols-1]) - u); 