% this benchmark specifies target rank we're going for
% computes the number of Krylov iterations based on the target rank
% computes the rank-target approximation
% block size is specified by user
function[] = RBKI_benchmark_Riley_analysis()
    %Mat1_RBKI_speed_comp_m_100000_n_100000_b_sz_start_16_b_sz_stop_128_num_matmuls_start_2_num_matmuls_stop_20
    Data_in = readmatrix("DATA_in/2024_03_22_runs/Mat1_RBKI_speed_comp_m_100000_n_100000_b_sz_start_16_b_sz_stop_128_num_matmuls_start_2_num_matmuls_stop_50.dat");
    
    for i = 1:size(Data_in, 1)
        % Mat 1 SVD: 1860675838
        % Mat 2 SVD: 1966542609
        % Mat 3 SVD: 1824193427
        % Mat 4 SVD: 2629361560
        % Mat 5 SVD: 2587508196
        Data_in(i, 8) = 1860675838;
    end
    
    Data_out = data_preprocessing(Data_in, 16, 128, 2, 50, 2);
    plot_speed_iters(Data_out, 16, 128, 2, 50);

    %file << b_sz << ",  " << RBKI.max_krylov_iters <<  ",  " << target_rank << ",  " << custom_rank << ",  " << residual_err_target <<  ",  " << residual_err_custom <<  ",  " << dur_rbki  << ",  " << dur_svd << ",\n";
end

function[Data_out] = data_preprocessing(Data_in, min_b_sz, max_b_sz, min_matmuls, max_matmuls, numiters)
    
    Data_out = [];
    num_b_sizes      = log2(max_b_sz) - log2(min_b_sz) + 1;
    num_matmuls = max_matmuls - min_matmuls + 1;
    
    i = 1;
  
    while i < num_b_sizes * num_matmuls * numiters
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

function[] = plot_speed_iters(Data, min_b_sz, max_b_sz, min_matmuls, max_matmuls)
    
    tiledlayout(1,2)

    num_b_sizes      = log2(max_b_sz) - log2(min_b_sz) + 1;
    num_matmuls = max_matmuls - min_matmuls + 1;

    % Plot error vs #iters
    nexttile
    for i = 1 : num_b_sizes
        % number of iters performed = (target_rank / b_sz) * 2
        y = Data((i - 1) * num_matmuls + 1 : i * num_matmuls, 6); % or 6 for custom
        x = Data((i - 1) * num_matmuls + 1 : i * num_matmuls, 2);
        
        plot(x, y, '-o', MarkerSize=10, LineWidth=2);
        set(gca,'YScale','log');
        hold on
        legend_entries{i} = ['B_{sz}=', num2str(Data(i * num_matmuls, 1))]; %#ok<AGROW>
    end
    grid on
    xlabel('#GEMM(A)', 'FontWeight','bold') 
    ylabel('sqrt(||AV - SigmaU||^2 + ||A^TU-VSigma||^2)', 'FontWeight','bold')
    legend(legend_entries);

    % Plot error vs speedup over SVD
    nexttile
    for i = 1 : num_b_sizes
        y = Data((i - 1) * num_matmuls + 1 : i * num_matmuls, 6); % or 6 for custom
        x = Data((i - 1) * num_matmuls + 1 : i * num_matmuls, 8) ./ Data((i - 1) * num_matmuls + 1 : i * num_matmuls, 7);
        
        plot(x, y, '-o', MarkerSize=10, LineWidth=2);
        set(gca,'YScale','log');
        hold on
    end
    
    
    grid on
    xlabel('Speedup over SVD', 'FontWeight','bold') 
end

