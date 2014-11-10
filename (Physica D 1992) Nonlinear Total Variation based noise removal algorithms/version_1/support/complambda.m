function lambda = complambda(u,u0,h,sigma)
% compute lambda according to equation (2.9c) in the paper
epsilon = 10^-3;
% epsilon = 1;
den  = sqrt(backdiffx(u).^2 + backdiffy(u).^2 + epsilon^2); % denominator
num1 = backdiffx(u0) .* backdiffx(u); % numerator
num2 = backdiffy(u0) .* backdiffy(u); % numerator

% Note that we must use mean(mean( )) here rather that sum(sum( ))
% The equation (2.9c) is a little misleading for us.
% However, if you look at equation (2.2c) and (2.3c), you would find that 
% the integration here is like computing the expectation.
lambda = -h/(2*sigma^2) * mean(mean(den - num1./den - num2./den));
