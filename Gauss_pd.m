% Gauss-LS (A not necessarily positive definite)
%theoretical part : B=A'*A, sigma = A*A'
function [x, log_resid, log_x, log_time, log_flops, n_iter] = Gauss_pd(A, b, e, max_time, verbose)
    [m, n] = size(A);

    [log_resid, log_x, log_time, log_flops] = deal([]);

    x = zeros(n, 1);
    n_iter = 0;
    time = 0;
    residue = norm(A * x - b) / norm(b);

    id = eye(n);
    mu = zeros(n,1);

    flops(0);
    tic;
    
    while residue > e && time < max_time
        n_iter = n_iter + 1;
        eta = mvnrnd(mu, id);
        eta = eta';
        x = x - dot(eta, A*x - b) / (dot(A*eta,eta)) * eta;
        residue = norm(A * x - b) / norm(b);

        if verbose == true
            if issparse(A) == true
                addflops(flops_spmul(A,x));
                addflops(flops_spmul(A,eta));
            else
                addflops(flops_mul(A,x));
                addflops(flops_mul(A,eta));
            end
            addflops(4*n); % for the dot products 
            % twice (n for multiplication, n for additions)
            addflops(flops_sqrt); %to compute the square root in normetaA
            addflops(n); % for the substraction in the dot product
            addflops(flops_div); % for the division 
            addflops(n); % for the scalar multiplication for myvec
            addflops(n); % for the last substraction with x
            log_resid = [log_resid, residue];
            log_x = [log_x, x];
            log_time = [log_time, toc];
            log_flops = [log_flops, flops];
        end
        time = toc;
    end
end
