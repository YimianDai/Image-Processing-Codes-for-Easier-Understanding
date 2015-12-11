function backDiffUy = backdiffy(u)
% backward diff in y direction in equation (2.9a)
% \Delta^y_{+}u_{ij} = u_{i,j+1} - u_{ij} 

cols = size(u,2);
% u([2:rows rows], :) is u_{i+1,j} in (2.9a)
backDiffUy = u(:, [2:cols cols]) - u; 