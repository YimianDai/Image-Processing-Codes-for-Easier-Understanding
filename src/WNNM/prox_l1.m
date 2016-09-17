function X_hat = prox_l1(X, tau)
    X_hat = max(abs(X) - tau, 0) .* sign(X);
end