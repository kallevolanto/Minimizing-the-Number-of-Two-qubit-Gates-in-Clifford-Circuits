function F = F_type_H(v)
    % applies a hadamard gate to each qubit index of v that is 1
    n = length(v);

    F = gf([eye(n)-diag(v),    diag(v);
            diag(v),      eye(n)-diag(v)]);
end