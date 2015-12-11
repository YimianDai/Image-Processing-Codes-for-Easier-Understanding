function backDiffUx = backdiffx(u)
% backward diff in x direction in equation (2.9a)
% \Delta^x_{+}u_{ij} = u_{i+1,j} - u_{ij} 

rows = size(u,1);
% u([2:rows rows], :) is u_{i+1,j} in (2.9a)
backDiffUx = u([2:rows rows], :) - u; 