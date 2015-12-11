function output = minmod(a,b)
% minmod is the minmod function in (2.9b) in the paper
% m(a,b) = minmod(a,b) = (\frac{sgn(a)+sgn(b)}{2}) min(|a|,|b|)

output = (sign(a) + sign(b)) / 2 * min(abs(a), abs(b));