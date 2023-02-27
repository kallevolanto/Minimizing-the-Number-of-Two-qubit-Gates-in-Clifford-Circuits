function N2 = bruhat_gate_decomposition(F, m)
    
    % Maslov, Dmitri, and Martin Roetteler. 
    % "Shorter stabilizer circuits via Bruhat decomposition 
    % and quantum circuit transformations." 
    % IEEE Transactions on Information Theory 64.7 (2018): 4729-4738.

    % Only the number of two-qubit gates is returned to save time
    % m is the parameter for the O(n^2/log(n)) algorithm used to decompose
    % the CNOT-type operators
    n = size(F,1)/2;
    omega = omega_mat(n);
    assert(isequal(F*omega*F',omega))
    
    if m > n
        % m must be smaller than n
        N2 = inf;
        return
    end

    % STEP 1
    % Obtain L, U, sigma and tau using Lemma 3
    C = F(n+1:2*n, 1:n);
    D = F(n+1:2*n, n+1:2*n);
    
    % Triangular matrices
    [L1,U1,P1,k] = lup_decomp(C);

    P1D1_mat = L1 *[C*U1 D*inv(U1)'];
   

    zero_cols = find(all(P1 == 0,1));
    zero_rows = find(all(P1 == 0,2));
    one_cols = find(~all(P1 == 0,1));
    one_rows = find(~all(P1 == 0,2));

    D1 = P1D1_mat(1:n, n+1:2*n);
    E = D1(zero_rows, zero_cols);
    [L2,L3,P2,n_minus_k] = llp_decomp(E);
    L2_prime = gf(zeros(n));
    L2_prime(one_rows, one_rows) = gf(eye(k));
    L2_prime(zero_rows, zero_rows) = L2;
    

    L3_prime = gf(zeros(n));
    L3_prime(one_cols, one_cols) = gf(eye(k));
    L3_prime(zero_cols, zero_cols) = L3;


    P1D2_mat = L2_prime * P1D1_mat * blkdiag(inv(L3_prime)', L3_prime);

    % Permutation matrices
    P1_k_times_k = P1(one_rows, one_cols);

    sigma = gf(zeros(n));
    sigma(1:k, one_rows) = P1_k_times_k';
    sigma(k+1:n, zero_rows) = gf(eye(n-k));
    
    % The ones are now ordered
    tau_1 = gf(zeros(n));
    tau_1(one_cols, 1:k) = gf(eye(k));
    tau_1(zero_cols, k+1:n) = gf(eye(n-k));

    % The left side is now an identity
    final_CD_mat = sigma * P1D2_mat * blkdiag(tau_1, tau_1);

    %
    P3 = final_CD_mat(k+1:n, n+k+1:2*n);
    tau_2 = blkdiag(gf(eye(k)), P3');

    L = L2_prime*L1;
    U = U1*inv(L3_prime)';
    %sigma = sigma
    tau = tau_1 * tau_2;



    % Part 2 
    % Complete the bruhat decomposition
    F2 = F_type1(sigma) * F_type1(inv(L')) * F * F_type1(U) * F_type1(tau);

    A1 = F2(1:k, 1:k);
    A3 = F2(k+1:n, 1:k);

    Q1 = [A1, A3';
          A3, zeros(n-k)];
    F3 = F_type2(Q1) * F2;

    D1 = F3(n+1:n+k, n+1:n+k);
    D2 = F3(n+1:n+k, n+k+1:2*n);
    B4_prime = F3(k+1:n, n+k+1:2*n);
    Q2 = [D1, D2;
          D2', B4_prime];
    F4 = F3 * F_type2(Q2);

    v = [ones(1,k), zeros(1,n-k)];
    
    % Check that the final result is a set of hadamard gates
    assert(isequal(F4, F_type_H(v)))
   
    % Set Q_prime1 and Q_prime2
    Q1_prime = sigma' * Q1 * sigma;
    Q2_prime = tau * Q2 * tau';
   
    % Check the decomposition is correct
    assert(isequal(F, F_type1(L') * F_type2(Q1_prime) * F_type1(sigma') * F_type_H(v) * F_type1(tau') * F_type2(Q2_prime) * F_type1(inv(U))))


    % STEP 3 
    % Obtain upper triangular U1 and U2 and diagonal Lambda1 and Lambda2
    % such that Q1_prime = U1*U1^T + Lambda1 and Q1_prime = U2*U2^T + Lambda2
    [U1, Lambda1] = upper_triangular_decomp(Q1_prime);
    [U2, Lambda2] = upper_triangular_decomp(Q2_prime);

    % STEP 4
    % Get the upper triangular cnot operators
    tri1 = L' * U1;
    tri2 = inv(U1);

    tri3 = U2;
    tri4 = inv(U2)*inv(U);

    % STEP 5
    % for the purpose of the thesis, only the amount of two-qubit gates are
    % needed, so only the triangular CNOT-type operators
    % defined by tri1, tri2, tri3 and tri4 need to be decomposed
    [N2_tri1,  ~] = CNOT_synth(tri1, m);
    [N2_tri2, ~] = CNOT_synth(tri2, m);
    [N2_tri3, ~] = CNOT_synth(tri3, m);
    [N2_tri4, ~] = CNOT_synth(tri4, m);

    N2 = N2_tri1 + N2_tri2 + N2_tri3 + N2_tri4;



end