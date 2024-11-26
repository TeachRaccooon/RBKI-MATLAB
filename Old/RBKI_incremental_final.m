% An algorithm that Rob proposed on 08/04/2023, supposed to be a rival to
% svdsketch().
function[U, Sigma, V, vecnorms_data_1, vecnorms_data_2, vecnorms_data_3, t_overhead] = RBKI_incremental_final(A, k, tol, maxiters, profiling, custom_rank)
    [m, n] = size(A);
    norm_A = norm(A, 'fro');
    sq_tol = tol^2;
    vecnorms_data_1 = [];
    vecnorms_data_2 = [];
    vecnorms_data_3 = [];
    t_overhead = 0;

    %Y_i = randn(n, k);
    Y_i = readmatrix("DATA_in/SKETCHING_OPERATOR.txt");
    Y_i = Y_i';
    size(Y_i)

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
            temp3 =  sqrt(temp1.^2 + temp2.^2);
    
            if size(vecnorms_data_1, 2) ~= size(temp1, 2)
                vecnorms_data_1 = [vecnorms_data_1, ones(size(vecnorms_data_1, 1), size(temp1, 2) - size(vecnorms_data_1, 2))];
                vecnorms_data_2 = [vecnorms_data_2, ones(size(vecnorms_data_2, 1), size(temp2, 2) - size(vecnorms_data_2, 2))];
                vecnorms_data_3 = [vecnorms_data_3, ones(size(vecnorms_data_3, 1), size(temp3, 2) - size(vecnorms_data_3, 2))];
            end
            vecnorms_data_1 = [vecnorms_data_1; temp1];
            vecnorms_data_2 = [vecnorms_data_2; temp2];
            vecnorms_data_3 = [vecnorms_data_3; temp2];
            t_overhead = t_overhead + toc(profiling_start);
        end

        i = i + 1;
    end

    fprintf("Total iters %d\n", i);
    writematrix(S, 'DATA_in/S_MATLAB.txt');
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
    
    U_RandLAPACK = readmatrix("DATA_in/U.txt")';
    VT_RandLAPACK = readmatrix("DATA_in/VT.txt")';
    S_RandLAPACK = readmatrix("DATA_in/S.txt");
    S_to_decomp_RandLAPACK = readmatrix("DATA_in/S_TO_DECOMPOSE.txt")';
    X_ev_RandLAPACK = readmatrix("DATA_in/X_ev.txt")';
    Y_od_RandLAPACK = readmatrix("DATA_in/Y_od.txt")';
    U_hat_RandLAPACK = readmatrix("DATA_in/U_hat.txt")';
    VT_hat_RandLAPACK = readmatrix("DATA_in/VT_hat.txt")';
   
    S_to_decomp_RandLAPACK = S_to_decomp_RandLAPACK(1:20, 1:16);
    fprintf("SVD_RandLAPACK %e\n", norm(S_to_decomp_RandLAPACK - U_hat_RandLAPACK * diag(S_RandLAPACK) * VT_hat_RandLAPACK));

    size(Sigma)
    size(S_RandLAPACK)
    fprintf("Sigma diff: %d\n", norm(Sigma - S_RandLAPACK, 'fro'));
    fprintf("U diff: %d\n", norm(U - U_RandLAPACK, 'fro'));
    fprintf("VT diff: %d\n", norm(V' - VT_RandLAPACK, 'fro'));
    fprintf("S diff: %d\n", norm(S - S_to_decomp_RandLAPACK, 'fro'));
    fprintf("U hat diff: %d\n", norm(U_hat - U_hat_RandLAPACK, 'fro'));
    fprintf("VT hat diff: %d\n", norm(V_hat' - VT_hat_RandLAPACK, 'fro'));
    fprintf("X diff: %d\n", norm(X_ev(:, 1:size(U_hat, 1)) - X_ev_RandLAPACK, 'fro'));
    fprintf("Y diff: %d\n", norm(Y_od(:, 1:size(V_hat, 1)) - Y_od_RandLAPACK(:, 1:size(V_hat, 1)), 'fro'));

   
    

    V_RandLAPACK = VT_RandLAPACK';
    S_RandLAPACK = diag(S_RandLAPACK);
    % Check A * V - U * Sigma
    residual_metric_target_1 = norm(A * V_RandLAPACK - U_RandLAPACK * S_RandLAPACK, 'fro'); 
    residual_metric_custom_1 = norm(A * V_RandLAPACK(:, 1:custom_rank) - U_RandLAPACK(:, 1:custom_rank) * S_RandLAPACK(1:custom_rank, 1:custom_rank), 'fro');

    % Check A' * U - V * Sigma
    residual_metric_target_2 = norm(A' * U_RandLAPACK - V_RandLAPACK * S_RandLAPACK, 'fro'); 
    residual_metric_custom_2 = norm(A' * U_RandLAPACK(:, 1:custom_rank) - V_RandLAPACK(:, 1:custom_rank) * S_RandLAPACK(1:custom_rank, 1:custom_rank), 'fro');
    size(A' * U_RandLAPACK(:, 1:custom_rank) - V_RandLAPACK(:, 1:custom_rank) * S_RandLAPACK(1:custom_rank, 1:custom_rank))


    fprintf("RandLAPACK\n");
    fprintf("Target residual %.20e\n", sqrt(residual_metric_target_1^2 + residual_metric_target_2^2));
    fprintf("Custom residual %.20e\n\n", sqrt(residual_metric_custom_1^2 + residual_metric_custom_2^2));
    %fprintf("Target residual %.20e\n", residual_metric_target_1);
    %fprintf("Target residual %.20e\n", residual_metric_target_2);
    %fprintf("Custom residual %.20e\n", residual_metric_custom_1);
    %fprintf("Custom residual %.20e\n\n", residual_metric_custom_2);
end