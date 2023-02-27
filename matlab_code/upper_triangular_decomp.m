function [U, Lambda] = upper_triangular_decomp(S)

    % S. Aaronson and D. Gottesman, “Improved simulation of stabilizer circuits,
    % ” Physical Review A, vol. 70, no. 5, p. 052328, 2004.

    % S is assumed to be symmetric, and corresponds to a symplectic matrix
    % F = [I S]
    %     [0 I]
    % Returns upper triangular U and diagonal Lambda
    % such that S = U U^T + Lambda

    n = size(S,1);
    S_copy = gf(S);
    U = gf(zeros(n));

    for i = n:-1:1
        
        v = S_copy(i, :);
        v(i) = 1;
        S_copy = S_copy + v'*v;

        % i:th column is v^T
        U(:,i) = v';

    end
    
    % Check that V'*V and S differ only by phase gates
    Lambda = U*U'+S;
    assert(isequal(diag(diag(Lambda)),Lambda))

    assert(isequal(triu(U), U))

end