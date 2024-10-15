% this benchmark specifies target rank we're going for
% computes the number of Krylov iterations based on the target rank
% computes the rank-target approximation
% block size is specified by user
function[] = RBKI_benchmark_Riley_analysis_target_rank()

    Data_in = readmatrix("DATA_in/2024_03_01_runs_target_rank/Mat1_old.txt");
    Data_out = data_preprocessing(Data_in, 8, 128, 256, 2048, 3);
    plot_speed_iters(Data_out, 8, 128, 256, 2048);

    %file << b_sz << ",  " << RBKI.max_krylov_iters <<  ",  " << target_rank << ",  " << custom_rank << ",  " << residual_err_target <<  ",  " << residual_err_custom <<  ",  " << dur_rbki  << ",  " << dur_svd << ",\n";
end

function[Data_out] = data_preprocessing(Data_in, min_b_sz, max_b_sz, min_target_rank, max_traget_rank, numiters)
    
    Data_out = [];
    num_b_sizes      = log2(max_b_sz) - log2(min_b_sz) + 1;
    num_target_ranks = log2(max_traget_rank) - log2(min_target_rank) + 1;
    
    i = 1;
  
    while i < num_b_sizes * num_target_ranks * numiters
        best_speed = intmax;
        best_speed_idx = 0;
        for j = 1:numiters
            if Data_in(i, 7) < best_speed
                best_speed = Data_in(i, 7);
                best_speed_idx = i;
            end
            i = i + 1;
        end
        Data_out = [Data_out; Data_in(best_speed_idx, :)]; %#ok<AGROW>
    end
end

function[] = plot_speed_iters(Data, min_b_sz, max_b_sz, min_target_rank, max_traget_rank)
    
    tiledlayout(1,2)

    num_b_sizes      = log2(max_b_sz) - log2(min_b_sz) + 1;
    num_target_ranks = log2(max_traget_rank) - log2(min_target_rank) + 1;

    % Plot error vs #iters
    nexttile
    for i = 1 : num_b_sizes
        % number of iters performed = (target_rank / b_sz) * 2
        y = Data((i - 1) * num_target_ranks + 1 : i * num_target_ranks, 5); % or 6 for custom
        x = Data((i - 1) * num_target_ranks + 1 : i * num_target_ranks, 2);
        
        plot(x, y, '-o', MarkerSize=10, LineWidth=2);
        set(gca,'YScale','log');
        hold on
        legend_entries{i} = ['B_{sz}=', num2str(Data(i * num_target_ranks, 1))]; %#ok<AGROW>
    end
    grid on
    xlabel('#GEMM(A)', 'FontWeight','bold') 
    ylabel('sqrt(||AV - SigmaU||^2 + ||A^TU-VSigma||^2)', 'FontWeight','bold')
    legend(legend_entries);

    % Plot error vs speedup over SVD
    nexttile
    for i = 1 : num_b_sizes
        y = Data((i - 1) * num_target_ranks + 1 : i * num_target_ranks, 5); % or 6 for custom
        x = Data((i - 1) * num_target_ranks + 1 : i * num_target_ranks, 8) ./ Data((i - 1) * num_target_ranks + 1 : i * num_target_ranks, 7);
        
        plot(x, y, '-o', MarkerSize=10, LineWidth=2);
        set(gca,'YScale','log');
        hold on
    end
    
    
    grid on
    xlabel('Speedup over SVD', 'FontWeight','bold') 
end