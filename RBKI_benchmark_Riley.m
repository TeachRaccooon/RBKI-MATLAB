% this benchmark specifies target rank we're going for
% computes the number of Krylov iterations based on the target rank
% computes the rank-target approximation
% block size is specified by user
function[] = RBKI_benchmark_Riley()

    Data_out = prepare_data()
    plot_3d(Data_out, 32, 128, 256, 1024, '||AV-SU||')
    hold off
    figure()
    plot_3d(Data_out, 32, 128, 256, 1024, '||A^TU-VS||')
    
end

function[] = plot_2d(Data, min_b_sz, max_b_sz, min_target_rank, max_traget_rank, accuracy_name)

    legend_entries = [];
    if strcmp(accuracy_name, '||AV-SU||')
        y_data_col = 3;
    elseif strcmp(accuracy_name, '||A^TU-VS||')
        y_data_col = 4;
    end

    num_b_sizes      = log2(max_b_sz) - log2(min_b_sz) + 1;
    num_target_ranks = log2(max_traget_rank) - log2(min_target_rank) + 1;

    for i = 1 : num_b_sizes
        x = Data((i - 1) * num_target_ranks + 1 : i * num_target_ranks, 5);
        y = Data((i - 1) * num_target_ranks + 1 : i * num_target_ranks, y_data_col);
        plot(x, y, '-o', MarkerSize=10, LineWidth=2);
        set(gca,'YScale','log');
        hold on
        legend_entries{i} = ['B_{sz}=', num2str(Data(i * num_target_ranks, 1))]; %#ok<AGROW>
    end
    
    
    grid on

    title('RBKI benchmark')
    xlabel('Speedup over SVD') 
    ylabel(accuracy_name)
    legend(legend_entries);
end

function[] = plot_3d(Data, min_b_sz, max_b_sz, min_target_rank, max_traget_rank, accuracy_name)

    legend_entries = [];
    if strcmp(accuracy_name, '||AV-SU||')
        y_data_col = 3;
    elseif strcmp(accuracy_name, '||A^TU-VS||')
        y_data_col = 4;
    end

    num_b_sizes      = log2(max_b_sz) - log2(min_b_sz) + 1;
    num_target_ranks = log2(max_traget_rank) - log2(min_target_rank) + 1;

    for i = 1 : num_b_sizes
        x = Data((i - 1) * num_target_ranks + 1 : i * num_target_ranks, 5);
        y = Data((i - 1) * num_target_ranks + 1 : i * num_target_ranks, y_data_col);
        z = Data((i - 1) * num_target_ranks + 1 : i * num_target_ranks, 2);
        
        plot3(x, y, z, '-o', MarkerSize=10, LineWidth=2);
        set(gca,'YScale','log');
        hold on
        legend_entries{i} = ['B_{sz}=', num2str(Data(i * num_target_ranks, 1))]; %#ok<AGROW>
    end
    
    
    grid on

    title('RBKI benchmark')
    xlabel('Speedup over SVD') 
    zlabel('Target rank') 
    ylabel(accuracy_name)
    legend(legend_entries);
end


function [Data_out] = call_RBKI(A, b_sz, tol, target_rank, Data_out, U_svd, Sigma_svd, V_svd, t_svd)

    numiters = (target_rank * 2) / b_sz;
    %numiters = ceil((size(A, 2) * 2) / b_sz); % This is the absolute maximum amount of iterations we can perform

    rbki_start = tic;
    [U_rbki, Sigma_rbki, V_rbki, vecnorms_data_1, vecnorms_data_2, t_overhead] = RBKI_incremental_final(A, b_sz, tol, numiters, 0);
    t_rbki = toc(rbki_start) - t_overhead;

    %hold off
    %figure()
    %plot_residuals(vecnorms_data_1, b_sz, '||AV-SU||');
    %hold off
    %figure()
    %plot_residuals(vecnorms_data_2, b_sz, '||A^TU-VS||');


    Sigma_rbki = diag(Sigma_rbki);
    target_rank = min(target_rank, size(Sigma_rbki, 1));
    custom_rank = ceil(target_rank / 8);

    % Check A * V - U * Sigma
    residual_metric_target_1 = norm(A * V_rbki - U_rbki * Sigma_rbki, 'fro') / sqrt(target_rank); 
    residual_metric_custom_1 = norm(A * V_rbki(:, 1:custom_rank) - U_rbki(:, 1:custom_rank) * Sigma_rbki(1:custom_rank, 1:custom_rank), 'fro') / sqrt(custom_rank);
    fprintf("Tartget residual  ||AV-SU||: %.20e\n", residual_metric_target_1);
    fprintf("Custom residual ||A^TU-VS||: %.20e\n", residual_metric_custom_1);


    % Check A' * U - V * Sigma
    residual_metric_target_2 = norm(A' * U_rbki - V_rbki * Sigma_rbki, 'fro') / sqrt(target_rank); 
    residual_metric_custom_2 = norm(A' * U_rbki(:, 1:custom_rank) - V_rbki(:, 1:custom_rank) * Sigma_rbki(1:custom_rank, 1:custom_rank), 'fro') / sqrt(custom_rank);
    fprintf("Tartget residual  ||AV-SU||: %.20e\n", residual_metric_target_2);
    fprintf("Custom residual ||A^TU-VS||: %.20e\n", residual_metric_custom_2);

    Data_out = [Data_out; b_sz, target_rank, residual_metric_target_1, residual_metric_target_2, t_svd / t_rbki];
end

function[] = plot_residuals(vecnorms_data, b_sz, name)

    % Plots norms of the columns of the resigual error matrix.
    % At every iteration, there are more columns to consider.
    size(vecnorms_data, 1)
    
    % This variable controls how much stuff we actually plot.
    % Anything after "data_stop" is garbage.
    data_stop = b_sz;
    % On the plot below for mat 1, odd iterations produce the "wave" plot.
    for i = 1:size(vecnorms_data, 1)
        if mod(i, 2) == 0
            semilogy(vecnorms_data(i, 1:data_stop),'Color', [1, 0, 0, 0.2]);
            data_stop = data_stop + b_sz;
        else 
            semilogy(vecnorms_data(i, 1:data_stop),'Color', [0, 0, 1, 0.2]);
        end
        hold on
    end
    % Legend & its transparency
    [~,BLicons] = legend('Odd Iter','Even Iter');
    PatchInLegend = findobj(BLicons, 'type', 'patch');
    set(PatchInLegend, 'facea', 1)
    % Title
    title(['Column norms of the residual error matrix (' name ')']);
end

function[Data_out] = prepare_data()
    A = readmatrix("DATA_in/test_mat_10k_rank_2k/RBKI_test_mat1.txt");
    [m, n] = size(A);
    tol = 2.5119e-14;

    b_sz = 32;
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
