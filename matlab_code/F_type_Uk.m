function F = F_type_Uk(n,k)
    % k has to be between 1 and n
    U_k = blkdiag(eye(k), zeros(n-k));
    L_n_k = blkdiag(zeros(k), eye(n-k));

    F = gf([U_k, L_n_k;
            L_n_k, U_k]);
end