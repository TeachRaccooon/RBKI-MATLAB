% this benchmark specifies target rank we're going for
% computes the number of Krylov iterations based on the target rank
% computes the rank-target approximation
% block size is specified by user
function[] = RBKI_benchmark_Riley()

    Data_out = prepare_data();
end


function [Data_out] = call_RBKI(A, b_sz, tol, target_rank, Data_out, U_svd, Sigma_svd, V_svd, t_svd)

    %numiters = (target_rank * 4) / b_sz;
    numiters = ceil((size(A, 2) * 2) / b_sz); % This is the absolute maximum amount of iterations we can perform

    tic;
    [U_rbki, Sigma_rbki, V_rbki, vecnorms_data_1, vecnorms_data_2] = RBKI_incremental_final(A, b_sz, tol, numiters);
    t_rbki = toc;

    
    % Plots norms of the columns of the resigual error matrix.
    % At every iteration, there are more columns to consider.
    size(vecnorms_data_1, 1)
    
    % This variable controls how much stuff we actually plot.
    % Anything after "data_stop" is garbage.
    data_stop = b_sz;
    % On the plot below for mat 1, odd iterations produce the "wave" plot.
    for i = 1:size(vecnorms_data_1, 1)
        if mod(i, 2) == 0
            semilogy(vecnorms_data_1(i, 1:data_stop),'Color', [1, 0, 0, 0.2]);
            data_stop = data_stop + b_sz;
        else 
            semilogy(vecnorms_data_1(i, 1:data_stop),'Color', [0, 0, 1, 0.2]);
        end
        hold on
    end
    % Legend & its transparency
    [BL,BLicons] = legend('Odd Iter','Even Iter')
    PatchInLegend = findobj(BLicons, 'type', 'patch');
    set(PatchInLegend, 'facea', 1)
    % Title
    title('Column norms of the residual error matrix (||AV-SU||)')


    hold off
    figure()

    data_stop = b_sz;
    % On the plot below for mat 1, odd iterations produce the "wave" plot.
    for i = 1:size(vecnorms_data_2, 1)
        if mod(i, 2) == 0
            semilogy(vecnorms_data_2(i, 1:data_stop),'Color', [1, 0, 0, 0.2]);
            data_stop = data_stop + b_sz;
        else 
            semilogy(vecnorms_data_2(i, 1:data_stop),'Color', [0, 0, 1, 0.2]);
        end
        hold on
    end
    % Legend & its transparency
    [BL,BLicons] = legend('Odd Iter','Even Iter')
    PatchInLegend = findobj(BLicons, 'type', 'patch');
    set(PatchInLegend, 'facea', 1)
    % Title
    title('Column norms of the residual error matrix (||A^TU-VS||)')


    Sigma_rbki = diag(Sigma_rbki);
    target_rank = min(target_rank, size(Sigma_rbki, 1));
    custom_rank = ceil(target_rank / 8);

    % Check A * V - U * Sigma
    residual_metric_target = norm(A * V_rbki - U_rbki * Sigma_rbki, 'fro') / sqrt(target_rank); 
    residual_metric_custom = norm(A * V_rbki(:, 1:custom_rank) - U_rbki(:, 1:custom_rank) * Sigma_rbki(1:custom_rank, 1:custom_rank), 'fro') / sqrt(target_rank);
    fprintf("Tartget residual: %.20e\n", residual_metric_target);
    fprintf("Custom residual: %.20e\n", residual_metric_custom);


%{
    VT_svd     = V_svd'; 
    VT_rbki    = V_rbki';
    Sigma_svd  = diag(Sigma_svd);
    Sigma_rbki = diag(Sigma_rbki);

    target_rank = min(target_rank, size(Sigma_rbki, 1))

    custom_rank = target_rank / 8;
    A_rbki_custom = U_rbki(:, 1:custom_rank) * Sigma_rbki(1:custom_rank, 1:custom_rank) * VT_rbki(1:custom_rank, :);
    A_svd_custom  = U_svd(:, 1:custom_rank) * Sigma_svd(1:custom_rank, 1:custom_rank) * VT_svd(1:custom_rank, :);

    A_rbki_target = U_rbki(:, 1:target_rank) * Sigma_rbki(1:target_rank, 1:target_rank) * VT_rbki(1:target_rank, :);
    A_svd_target  = U_svd(:, 1:target_rank) * Sigma_svd(1:target_rank, 1:target_rank) * VT_svd(1:target_rank, :);

    err_svals_rbki_target =  norm(Sigma_svd(1:target_rank, 1:target_rank)  - Sigma_rbki(1:target_rank, 1:target_rank), "fro") / norm(Sigma_svd(1:target_rank, 1:target_rank),  "fro");
    err_fro_rbki_target   =  norm(A_svd_target  - A_rbki_target, "fro") / norm(A_svd_target,  "fro");

    err_svals_rbki_custom =  norm(Sigma_svd(1:custom_rank, 1:custom_rank)  - Sigma_rbki(1:custom_rank, 1:custom_rank), "fro") / norm(Sigma_svd(1:custom_rank, 1:custom_rank),  "fro");
    err_fro_rbki_custom   =  norm(A_svd_custom  - A_rbki_custom, "fro") / norm(A_svd_custom,  "fro");

    % Check A * V - U * Sigma
    residual_metric_target = norm(A * V_rbki - U_rbki * Sigma_rbki, 'fro') / sqrt(target_rank); 
    residual_metric_custom = norm(A * V_rbki(:, 1:custom_rank) - U_rbki(:, 1:custom_rank) * Sigma_rbki(1:custom_rank, 1:custom_rank), 'fro') / sqrt(target_rank); 

    fprintf("Target ||S_svd  - S_rbki||_F/||S_svd||_F: %.20e\n", err_svals_rbki_target);
    fprintf("Target ||A_svd  - A_rbki||_F/||A_svd||_F: %.20e\n", err_fro_rbki_target);
    fprintf("Custom ||S_svd  - S_rbki||_F/||S_svd||_F: %.20e\n", err_svals_rbki_custom);
    fprintf("Custom ||A_svd  - A_rbki||_F/||A_svd||_F: %.20e\n", err_fro_rbki_custom);

    fprintf("Tartget residual: %.20e\n", residual_metric_target);
    fprintf("Custom residual: %.20e\n", residual_metric_custom);

    fprintf("RBKI is %ex faster than SVD\n\n", t_svd/ t_rbki);

    s_svd = diag(Sigma_svd(1:target_rank, 1:target_rank));
    s_rbki = diag(Sigma_rbki(1:target_rank, 1:target_rank));

    x = 1:target_rank;
    hold on
    semilogy(x, diag(Sigma_svd(1:target_rank, 1:target_rank)), 'Color',  'black','LineWidth', 1.6)
    hold on
    semilogy(x, diag(Sigma_rbki(1:target_rank, 1:target_rank)), 'Color',  'blue','LineWidth', 1.6)
    %hold off
    %semilogy(x, (abs(s_svd - s_rbki) + 10^-16) / norm(s_svd, 2))
    %legend('SVD','RBKI')

    %Data_out = [Data_out; b_sz, target_rank, err_svals_rbki, t_rbki, t_svd];
%}
end

function[Data_out] = prepare_data()
    A = readmatrix("DATA_in/test_mat_1k_rank_200/RBKI_test_mat2.txt");
    [m, n] = size(A);
    tol = 2.5119e-14;

    b_sz = 20;
    b_sz_max = 20;
    target_rank = 1000;
    target_rank_max = 1000;
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
