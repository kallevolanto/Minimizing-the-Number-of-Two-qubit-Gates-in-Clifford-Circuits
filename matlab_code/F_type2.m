function F = F_type2(S)
    % S has to be a symmetric n X n matrix
    S = gf(S);
    assert(isequal(S,S'))
    n = size(S,1);

    F = [eye(n), S;
        zeros(n), eye(n)];
end