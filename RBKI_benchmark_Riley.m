% this benchmark specifies target rank we're going for
% computes the number of Krylov iterations based on the target rank
% computes the rank-target approximation
% block size is specified by user
function[] = RBKI_benchmark_Riley()

    Data_out = prepare_data()
end


function [Data_out] = call_RBKI(A, b_sz, tol, target_rank, Data_out, U_svd, Sigma_svd, V_svd, t_svd)

    numiters = (target_rank * 2) / b_sz;

    tic;
    [U_rbki, Sigma_rbki, V_rbki] = RBKI_incremental_final(A, b_sz, tol, numiters);
    t_rbki = toc;

    VT_svd     = V_svd'; 
    VT_rbki    = V_rbki';
    Sigma_svd  = diag(Sigma_svd);
    Sigma_rbki = diag(Sigma_rbki);

    custom_rank = target_rank / 2;
    A_rbki_custom = U_rbki(:, 1:custom_rank) * Sigma_rbki(1:custom_rank, 1:custom_rank) * VT_rbki(1:custom_rank, :);
    A_svd_custom  = U_svd(:, 1:custom_rank) * Sigma_svd(1:custom_rank, 1:custom_rank) * VT_svd(1:custom_rank, :);

    A_rbki_target = U_rbki(:, 1:target_rank) * Sigma_rbki(1:target_rank, 1:target_rank) * VT_rbki(1:target_rank, :);
    A_svd_target  = U_svd(:, 1:target_rank) * Sigma_svd(1:target_rank, 1:target_rank) * VT_svd(1:target_rank, :);

    err_svals_rbki_target =  norm(Sigma_svd(1:target_rank, 1)  - Sigma_rbki(1:target_rank, 1), "fro") / norm(Sigma_svd(1:target_rank, 1),  "fro");
    err_fro_rbki_target   =  norm(A_svd_target  - A_rbki_target, "fro") / norm(A_svd_target,  "fro");

    err_svals_rbki_custom =  norm(Sigma_svd(1:custom_rank, 1)  - Sigma_rbki(1:custom_rank, 1), "fro") / norm(Sigma_svd(1:custom_rank, 1),  "fro");
    err_fro_rbki_custom   =  norm(A_svd_custom  - A_rbki_custom, "fro") / norm(A_svd_custom,  "fro");

    fprintf("Target ||S_svd  - S_rbki||_F/||S_svd||_F: %.20e\n", err_svals_rbki_target);
    fprintf("Target ||A_svd  - A_rbki||_F/||A_svd||_F: %.20e\n", err_fro_rbki_target);
    fprintf("Custom ||S_svd  - S_rbki||_F/||S_svd||_F: %.20e\n", err_svals_rbki_custom);
    fprintf("Custom ||A_svd  - A_rbki||_F/||A_svd||_F: %.20e\n", err_fro_rbki_custom);
    fprintf("RBKI is %ex faster than SVD\n\n", t_svd/ t_rbki);

    %Data_out = [Data_out; b_sz, target_rank, err_svals_rbki, t_rbki, t_svd];

end

function[Data_out] = prepare_data()
    A = readmatrix("DATA_in/test_mat_10k_rank_2k/RBKI_test_mat1.txt");
    [m, n] = size(A);
    tol = 2.5119e-14;

    b_sz = 8;
    b_sz_max = 128;
    target_rank = 256;
    target_rank_max = 1024;
    target_rank_start = target_rank;

    Data_out = [];

    tic;
    [U_svd, Sigma_svd, V_svd] = svd(A, 'econ', 'vector');
    t_svd = toc;

    while b_sz <= b_sz_max
        while target_rank <= target_rank_max
            Data_out = call_RBKI(A, b_sz, tol, target_rank, Data_out, U_svd, Sigma_svd, V_svd, t_svd);
            target_rank = target_rank * 2;
        end
        target_rank = target_rank_start;
        b_sz = b_sz * 2;
    end
end
