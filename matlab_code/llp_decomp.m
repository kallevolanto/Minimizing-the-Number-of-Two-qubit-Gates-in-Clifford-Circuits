function [L_left,L_right,P,k] = llp_decomp(M)
    
    % G. Strang, "Fast transforms: Banded matrices with banded inverses,” Proceedings
    % of the National Academy of Sciences, vol. 107, no. 28, pp. 12413–12416, 2010.

    % Return lower triangular matrices L_left and L_right
    % such that L_left*M*L_right = P
    % where P is a permutation pattern
    % i.e only has one non-zero entry per row or column
    
    n = size(M,1);
    assert(size(M,1)==size(M,2))
    M = gf(M);

    
    L_left = gf(eye(n));
    L_right = gf(eye(n));

    % Loop through each row
    for row = 1:n 
        curr_row = M(row, :);
        
        % Pick the rightmost non-zero element
        rightmost_col = find(curr_row==1,1,'last');
        if isempty(rightmost_col)
            % This row can be skipped as it is already zero
            continue
        end
        
        % Remove the non-zero elements from the right (column operations)
        L_right_update = gf(eye(n));
        L_right_update(rightmost_col, :) = curr_row;
        L_right = L_right * L_right_update;
        M = M * L_right_update;

        % Remove the non-zero elements elements below (row operations)
        L_left_update = gf(eye(n));
        L_left_update(:, row) = M(:, rightmost_col);
        L_left = L_left_update * L_left;
        M = L_left_update * M;
    end
    
    % The end result is the permutation pattern P 
    P = M;
    k = rank(P);
    assert(sum(P==1,'all')==k)
end