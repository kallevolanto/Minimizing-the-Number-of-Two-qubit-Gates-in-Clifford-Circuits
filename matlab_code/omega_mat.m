function omega = omega_mat(n)

    omega = gf([zeros(n), eye(n);
            eye(n), zeros(n)]);

end