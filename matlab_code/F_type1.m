function F = F_type1(P)
    % P has to be invertible n X n matrix
    P = gf(P);
    n = size(P,1);

    F = [P, zeros(n);
         zeros(n), inv(P')];
end