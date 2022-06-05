function [x, log_resid, log_x, log_time, log_flops, n_iter] = CD_LS(A, b, e, max_time, verbose)
    [m, n] = size(A);

    [log_resid, log_x, log_time, log_flops] = deal([]);

    x = zeros(n, 1);
    n_iter = 0;
    time = 0;
    residue = norm(A * x - b) / norm(b);
    I = eye(n); 

    frob_norm = norm(A,"fro"); 
    prob_array = vecnorm(A).^2 ./ frob_norm^2;
    if issparse(A) == true
        prob_array = full(prob_array);
    end
    dpdf = makedist('Multinomial', 'Probabilities', prob_array);

    flops(0);
    tic;
    
    while residue > e && time < max_time
        n_iter = n_iter + 1;
        n_col = random(dpdf);
        x = x - dot(A(:,n_col), (A*x - b) / (norm(A(:,n_col))^2)) .* I(:, n_col);

        residue = norm(A * x - b) / norm(b);

        if verbose == true
            if issparse(A) == true
                addflops(flops_spmul(A, x));
            else
                addflops(flops_mul(A, x));
            end

            addflops(n); % for the dot product
            addflops(n); % for the substraction
            addflops(n*flops_div); % for the division in each term
            addflops(2*n+2);  % to do this operation (norm(A(i,:))^2)
            addflops(n); % for the scalar multiplication for myvec

            log_resid = [log_resid, residue];
            log_x = [log_x, x];
            log_time = [log_time, toc];
            log_flops = [log_flops, flops];
        end
        time = toc;
    end
end
