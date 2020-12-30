function s = S(u)
    % substitution to handle the nonlinear terms
    alpha = 4;
    s = -u.*(1-u)+(u./(u+alpha));
end