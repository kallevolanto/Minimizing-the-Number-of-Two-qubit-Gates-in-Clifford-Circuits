function [N2, CNOT_list] = CNOT_synth(A, m)

    % Patel, Ketan N., Igor L. Markov, and John P. Hayes. 
    % "Optimal synthesis of linear reversible circuits." 
    % Quantum Inf. Comput. 8.3 (2008): 282-294.
    
    % N2 is the number of CNOTs
    % CNOT_list is a list of the CNOT pairs (control, target)
    
    % A needs to be invertible and triangular
    A = A==1;
    n = size(A,1);
    assert(gfrank(A,2)==n)
    assert(isequal(tril(A),A) || isequal(triu(A),A))
    
    % To skip the case where A is just 1
    if n==1
        N2 = 0;
        CNOT_list =[];
        return
    end

    upper = false;
    if isequal(triu(A),A)
        % make the matrix lower triangular
        upper = true;
        A = A';
    end

    % m is the partition size
    assert(m <=  n);
    
    CNOT_list = [];
    N2 = 0;

    % Add the target row to control row
    function update_A(control, target)
        A(control, :) = mod(A(control, :) + A(target, :),2);
        CNOT_list(end+1, :) = [target, control];
        N2 = N2 + 1;

        assert(control > target)
    end


    % Loop for each section (move diagonally)
    for section = 1:ceil(n/m)
        % Remove duplicate sub-rows
        % 0 means pattern is not yet found
        col_start = (m*(section-1)+1);
        col_end = (m*(section));
        if col_end > n
            col_end = n;
            patterns = zeros(2^(col_end-col_start+1), 1);
        else
            patterns = zeros(2^(m), 1);
        end
       
        % Loop through each row starting at the diagonal element
        row_start = (m*(section-1)+1);
        row_end = n;
        for row = row_start:row_end
            patt = A(row, col_start:col_end);
            patt_num = bi2de(patt)+1; % +1 cause matlab
            
            patt_row = patterns(patt_num);
            if patt_row == 0
                % 0 means pattern is not yet found
                patterns(patt_num) = row;
            else
                % duplicate pattern found, remove it
                update_A(row, patt_row)
            end
    
        end
        
        % Remove the remaining ones in this section
        for col = col_start:col_end
            
            % Check that there is a 1 on diagonal
            diag_one = true;
            if A(col, col) == 0
                diag_one = false;
            end
            
            for row = col+1:n
                if A(row, col) == 1
                    % If there is no 1 in the diagonal, put it there
                    % This never happens if A is lower triangular to begin
                    % with
                    if diag_one == false
                        update_A(col, row)
                        diag_one = true;
                    end

                    % Remove the 1 from row
                    update_A(row,col);

                end
            end
    

        end

    end


    assert(isequal(A, eye(n)))

    % If the matrix was upper triangular
    % The CNOTs have to be flipped and the order has to be reversed
    if upper
        flip(flip(CNOT_list,2),1);
    end


end