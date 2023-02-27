function F = rand_symp_mat(n)

    % T. Can, “An algorithm to generate a unitary transformation from logarithmically
    % many random bits," Research Independent Study, 2017.


    % Generates a uniformly random symplectic matrix
    
    % F = [A B
    %      C D]
    % [A B] is a random maximal symplectically orthogonal basis
    V = orthogonal_symp_basis(n);
    A = V(:, 1:n);
    B = V(:, n+(1:n));


    % Gaussian elimination to find M and N
    [~, M_orig, N, rank_A] = gf2rref(A);
    N = gf(N);
    M_orig = gf(M_orig);
    
    new_form = M_orig*[A, B]*F_type1(N);

    %  [A B] should now be
    %   [I_k    0   Q_k R ]
    %   [0      0   0   B_(n-k)]

    % Change M a little to make the result into
    %   [I_k    0   Q_k 0 ]
    %   [0      0   0   I_(n-k)]
    
    % (k = rank_A)
    
    % make the lower right part identity 
    B_rest = new_form(rank_A+1:end, n+rank_A+1:end);
    [~, M_up1, ~,  rank_rest] = gf2rref(B_rest);
    M_up1 = gf(blkdiag(eye(rank_A), M_up1));  
    assert(rank_A+rank_rest==n) % rank_A should be the n-rank_rest
    
    
    % Remove the ones from the right of Q_k
    B_up_right = new_form(1:rank_A, n+rank_A+1:end);
    M_up2 = gf([eye(rank_A), B_up_right;
                zeros(rank_rest, rank_A), eye(rank_rest)]);
    
    % Update M
    M = M_up2*M_up1*M_orig;
    
    % Remove Q_k to get
    %   [I_k    0   0 0 ]
    %   [0      0   0   I_(n-k)]
    Q_k = new_form(1:rank_A, n+(1:rank_A));
    Q = blkdiag(Q_k, zeros(rank_rest));
    
    
    % A final change to get
    %   [I_n 0]
    final_form = M*[A, B]*F_type1(N)*F_type2(Q)*F_type_Uk(n,rank_A);
    
    assert(isequal(final_form, gf([eye(n), zeros(n)])))
    

    % We can now generate any F with [A B] as the first n rows just by
    % generating a random symmetric matrix
    
    P = gf(randi([0,1], n, n));
    P = P+P'+diag(diag(P));
    k = rank_A;
    

    F = F_type1(inv(M))*F_type3(P)*F_type_Uk(n,k)*F_type2(Q)*F_type1(inv(N));
    assert(isequal(F(1:n, :), [A B]))
    
end