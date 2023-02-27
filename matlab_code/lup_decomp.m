function [L,U,P,k] = lup_decomp(M)
    
    % G. Strang, “Fast transforms: Banded matrices with banded inverses,” Proceedings
    % of the National Academy of Sciences, vol. 107, no. 28, pp. 12413–12416, 2010.

    % Return lower triangular L and upper triangular U
    % such that L*M*U = P
    % where P is a permutation pattern
    % i.e only has one non-zero entry per row or column

    n = size(M,1);
    assert(size(M,1)==size(M,2))
    M = gf(M);

    
    L = gf(eye(n));
    U = gf(eye(n));

    % Loop through each row
    for row = 1:n 
        curr_row = M(row, :);
        
        % Pick the leftmost non-zero element
        leftmost_col = find(curr_row==1,1);
        if isempty(leftmost_col)
            % This row can be skipped as it is already zero
            continue
        end
        
        % Remove the non-zero elements from the right (column operations)
        U_update = gf(eye(n));
        U_update(leftmost_col, :) = curr_row;
        U = U * U_update;
        M = M * U_update;

        % Remove the non-zero elements elements below (row operations)
        L_update = gf(eye(n));
        L_update(:, row) = M(:, leftmost_col);
        L = L_update * L;
        M = L_update * M;
    end
    
    % The end result is the permutation pattern P 
    P = M;
    k = rank(P);
    assert(sum(P==1,'all')==k)
end



