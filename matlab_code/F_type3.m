function F = F_type3(S)
    % S has to be a symmetric n X n matrix
    S = gf(S); 
    assert(isequal(S,S'))
    n = size(S,1);

    F = [eye(n), zeros(n);
        S, eye(n)];
end