function V = orthogonal_symp_basis(n)

    % Generates a random maximal orthogoanlly symplectic basis
    
    % Algorithm taken from
    % TRUNG CAN
    % ADVISOR: ROBERT CALDERBANK
    % SPRING 2017, DUKE UNIVERSITY

    % RESEARCH INDEPENDENT STUDY
    % AN ALGORITHM TO GENERATE A UNITARY TRANSFORMATION
    % FROM LOGARITHMICALLY MANY RANDOM BITS


    % the basis
    V = gf([]);
    
    % The space that is orthogonal to V
    U = gf([eye(2*n)]);
    
    rank_V = 0;
    omega = omega_mat(n);

    for i = 1:n
        % Generate v_i as a random combination until it doesn't belong to span(V)
        v_i = randi([0,1], [1,size(U,1)]) * U;
        while rank([V; v_i])~=rank_V+1
    
            v_i = randi([0,1], [1,size(U,1)]) * U;
        end
        V = [V; v_i];
        rank_V = rank_V+1;
     
       
        % Compute the new orthogonal space
        inner_prod = U*omega*v_i';
    
        % Already orthogonal
        P_i = U(inner_prod==0, :);
    
        % Needs some work
        Q_i = U(inner_prod==1, :);
    
        % Choose the first q, remove it and add to others
        % Q_i is now orthogonal to V
        len_Q = size(Q_i,1);
        q = Q_i(1, :);
        Q_updated = Q_i(2:end, :) + repmat(q, len_Q-1, 1);
    
        % U is just the combination of P_i and Q_updated
        U = [P_i; Q_updated];
    end
    
    % Assert that V is in fact self-orthogonal
    assert(all(V*omega*V','all')==0)
end