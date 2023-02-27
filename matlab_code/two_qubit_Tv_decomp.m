function [N2, N1, NSwaps] = two_qubit_Tv_decomp(F)
    
    % Takes in a symplectic binary matrix F
    % 
    % Decomposes F into swap gates, two-qubit transvections and one-qubit
    % transvections
    %
    % Note that this function only returns the amount of operators used
   
    N2 = 0; % Number of 2-qubit transvections
    N1 = 0; % Number of 1-qubit transvections
    NSwaps = 0; % Number of swap gates

    F = gf(F);
    n = size(F,1)/2;
    omega = omega_mat(n);
    I = gf(eye(2*n));
    
    assert(isequal(F'*omega*F,omega));
    
    
    % Gives the defined submatrix R_F(i,j)
    function R_mat = sub_mat(i,j)
        R_mat = [F(i,j), F(i, n+j); F(i+n,j), F(i+n, n+j)];
    end

    
    % Updating F -> T_v*F
    function apply_transvection(v)
        F = F + (omega*v')*(v*F);
        num_qubits = sum(v(1:n)==1|v(n+1:end)==1);
        assert(num_qubits > 0)
        assert(num_qubits < 3)
        if num_qubits==2
            N2 = N2+1;
        else
            N1 = N1+1;
        end
    end

    % Updating F -> SWAP(i,j)*F
    % swap gate swaps rows (i,n+i) between (j, n+j)
    function apply_swap(i,j)
        assert(i > 0)
        assert(i <= n)
        assert(j > 0)
        assert(j <= n)

        perm_order = 1:n;
        perm_order(i) = j;
        perm_order(j) = i;
        
        % F -> SWAP(i,j)*F
        F = F([perm_order, n+perm_order], :);

        NSwaps = NSwaps + 1;
    end

    % outer loop
    for i = 1:n
        % List all the invertible submatrices R(j,i), j>i
        inv_list = [];
        for j = i:n
            if rank(sub_mat(j,i))==2
                inv_list = [inv_list, j];
            end
        end
        % There should be an odd number of items
        %assert(mod(length(inv_list),2)==1)
        
        % STEP 1
        % make R(i,i) invertible
        if inv_list(1) ~= i
            j = inv_list(1);
            apply_swap(i,j)
            
            % Storing the swap operation (not done here)
            %sequence = ['SWAP_ij' sequence];
        end
        
        
        % The list of transvections to be stored in the sequence
        % Each row corresponds to a transvection
        %V = [];

        % STEP 2 
        % Loop through the remaining invertible pairs and make them
        % non-invertible
        for idx = 2:2:length(inv_list)
            j = inv_list(idx);
            k = inv_list(idx+1);

            R_ji = sub_mat(j,i);
            R_ki = sub_mat(k,i);

            % Choose a,b,c,d as given in the algorithm
            c = F(n+k,i);
            d = F(k,i);

            % [a,b] are chose so that [a b]R_ji + [c d]R_ki = [1 0]
            vec = ([1 0] + [c d]*R_ki) * inv(R_ji);
            a = vec(1);
            b = vec(2);
            
            % Apply the transvection as given in the algorithm
            v = gf(zeros(1,2*n));
            v(j) = a;
            v(j+n) = b;
            v(k) = c;
            v(k+n) = d;
            apply_transvection(v)
            
            % Store the transvection (NOT DONE HERE)
            %V = [v; V];

            %R_ji and R_ki are no longer invertible
            %assert(rank(sub_mat(j,i))==1)
            %assert(rank(sub_mat(k,i))==1)
        end
            
        % STEP 3
        % Loop through the submatrices R(k,i), k > i and make them 0
        for k = i+1:n
            
            R_ii = sub_mat(i,i);
            R_ki = sub_mat(k,i);
    
            if rank(R_ki)==1
                % choose a,b,c,d as given in the algorithm
                % get non-zero columns of R_ki
                non_zero_cols = [0 0];
                
                if (R_ki(1,1) == 1 || R_ki(2,1)==1)
                    non_zero_cols(1) = 1;
                    c = R_ki(2,1);
                    d = R_ki(1,1);
                end
                if (R_ki(1,2) == 1 || R_ki(2,2)==1)
                    non_zero_cols(2) = 1;
                    c = R_ki(2,2);
                    d = R_ki(1,2);
                end

                % [a b] R_ii = non_zero_cols
                vec = non_zero_cols * inv(R_ii);
                a = vec(1);
                b = vec(2);
                
                % Apply the transvection as given in the algorithm
                v = gf(zeros(1,2*n));
                v(i) = a;
                v(i+n) = b;
                v(k) = c;
                v(k+n) = d;
                apply_transvection(v)
                
                % Store the transvection (NOT DONE HERE)
                %sequence = [T_v; sequence]; 
            end

            %R_ki is now 0
            %assert(rank(sub_mat(k,i))==0)
        end

        % Step 4
        % make it so that R(i,i) = I_2 using Gaussian elimination on the 
        % 2 x 2 matrix
        if F(i,i)==0
            % swap rows i, n+1
            v = gf(zeros(1,2*n));
            v(i) = 1;
            v(n+i) = 1;
            apply_transvection(v)
        end
        % now F(i,i) = 1
        if F(n+i, i)==1
            % add row i to row n+i
            v = gf(zeros(1,2*n));
            v(i) = 1;
            apply_transvection(v)
        end
        % now either
        % R(i,i) = [1 1
        %           0 1]
        % or
        % R(i,i) = [1 0
        %           0 1]
        if F(i, n+i)==1
            % add row n+1 to row i
            v = gf(zeros(1,2*n));
            v(n+i) = 1;
            apply_transvection(v)
        end
    end


    % Final test to see that the decomposition is complete
    %F
    assert(isequal(F,I))

   
end