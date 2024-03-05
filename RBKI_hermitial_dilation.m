% An algorithm that Riley proposed on 02/21/2024, supposed to solve the residual error metric issue of RBKI.
function[U, Sigma, V, vecnorms_data_1, vecnorms_data_2, t_overhead] = RBKI_hermitial_dilation(A_base, k, tol, maxiters, profiling)
    
    A = [zeros(size(A_base)), A_base; A_base', zeros(size(A_base))];

    [m, n] = size(A);
    norm_A = norm(A, 'fro');
    sq_tol = tol^2;
    vecnorms_data_1 = [];
    vecnorms_data_2 = [];
    t_overhead = 0;

    

    Y_i = randn(n, k);


    [X_i, ~] = qr(A * Y_i, 0);
    R = []; S = [];
    X_ev = X_i; Y_od = zeros(n, 0);
    i = 1;
    while 1
        if mod(i, 2) ~= 0
            Y_i = A' * X_i;
            R_i = Y_od' * Y_i;
            Y_i = Y_i - Y_od * R_i;
            % Reorthog
            Y_i = Y_i - Y_od *  Y_od' * Y_i;
            [Y_i, R_ii] = qr(Y_i, 0);
            R = [R, R_i; zeros(size(R_ii, 1), size(R, 2)), R_ii]; %#ok<AGROW>
            Y_od = [Y_od, Y_i]; %#ok<AGROW>
            % Early termination - vectors exhausted
            if (abs(R(end)) <= sqrt(eps('double')))
                disp("Early termination 1");
                break;
            end
        else
            X_i = A * Y_i;
            S_i = X_ev' * X_i;
            X_i = X_i - X_ev * S_i;
            % Reorthog
            X_i = X_i - X_ev * X_ev' * X_i;
            [X_i, S_ii] = qr(X_i, 0);
            S = [S, S_i; zeros(size(S_ii, 1), size(S, 2)), S_ii]; %#ok<AGROW>
            X_ev = [X_ev, X_i]; %#ok<AGROW>
            % Early termination - vectors exhausted
            if (abs(S(end)) <= sqrt(eps('double')))
                disp("Early termination 2");
                break;
            end
        end

        % ||A - \hat(A)||_F < eps * ||A||_F
        if norm(R, 'fro') > sqrt(1 - sq_tol) * norm_A
            disp("TERMINATION 3");
            break;
        end
        if i >= maxiters
            break;
        end

        if profiling
            profiling_start = tic;
            if mod(i, 2) ~= 0
                [U_hat, Sigma, V_hat] = svd(R', 'econ', 'vector');
                U = X_ev(:, 1:size(U_hat, 1)) * U_hat;
                V = Y_od(:, 1:size(V_hat, 1)) * V_hat;
            else
                [U_hat, Sigma, V_hat] = svd(S, 'econ', 'vector');
                U = X_ev(:, 1:size(U_hat, 1)) * U_hat;
                V = Y_od(:, 1:size(V_hat, 1)) * V_hat;
            end
            
            temp1 =  vecnorm(A * V - U * diag(Sigma));
            temp2 =  vecnorm(A' * U -  V * diag(Sigma));
    
            if size(vecnorms_data_1, 2) ~= size(temp1, 2)
                vecnorms_data_1 = [vecnorms_data_1, ones(size(vecnorms_data_1, 1), size(temp1, 2) - size(vecnorms_data_1, 2))];
                vecnorms_data_2 = [vecnorms_data_2, ones(size(vecnorms_data_2, 1), size(temp2, 2) - size(vecnorms_data_2, 2))];
            end
            vecnorms_data_1 = [vecnorms_data_1; temp1];
            vecnorms_data_2 = [vecnorms_data_2; temp2];
            t_overhead = t_overhead + toc(profiling_start);
        end

        i = i + 1;
    end

    %fprintf("Total iters %d\n", i);

    size(R)
    size(S)
 
    if mod(i, 2) ~= 0
        %fprintf("SVD on R;\n")
        [U_hat, Sigma, V_hat] = svd(R', 'econ', 'vector');
        U = X_ev(:, 1:size(U_hat, 1)) * U_hat;
        V = Y_od(:, 1:size(V_hat, 1)) * V_hat;
    else
        %fprintf("SVD on S\n")
        [U_hat, Sigma, V_hat] = svd(S, 'econ', 'vector');
        U = X_ev(:, 1:size(U_hat, 1)) * U_hat;
        V = Y_od(:, 1:size(V_hat, 1)) * V_hat;
    end
    
    temp1 =  vecnorm(A * V - U * diag(Sigma))
    temp2 =  vecnorm(A' * U -  V * diag(Sigma))

    semilogy(temp1, 'Color', [1, 0, 0, 0.9])
    hold on
    semilogy(temp2, 'Color', [0, 1, 0, 0.9])
end