% This code is temporary bc idk how to properly incorporate the fact that
% some of the data for block sizes is missing
function[] = RBKI_benchmark_Riley_analysis_temporary()
    Data_in = readmatrix("DATA_in/2024_04_01_runs_large_matmuls/Mat5_RBKI_speed_comp_m_100000_n_100000_b_sz_start_1_b_sz_stop_128_num_matmuls_start_2_num_matmuls_stop_50.txt");
    
    for i = 1:size(Data_in, 1)
        % Mat 1 SVD: 1860675838
        % Mat 2 SVD: 1966542609
        % Mat 3 SVD: 1824193427
        % Mat 4 SVD: 2629361560
        % Mat 5 SVD: 2587508196
        Data_in(i, 8) = 1860675838;
    end
    
    max_b_sz = 128;
    min_b_sz = 1;

    max_matmuls = 50;
    min_matmuls = 2;

    % We don't have b_sz 4 in this dataset, so need to add -1.
    num_b_sizes = (log2(max_b_sz) - log2(min_b_sz) + 1) -1;
    num_matmuls = max_matmuls - min_matmuls + 1;


    Data_out = data_preprocessing(Data_in, num_b_sizes, num_matmuls, 2);
    plot_speed_iters(Data_out, num_b_sizes, num_matmuls, num_matmuls, 1);

    %file << b_sz << ",  " << RBKI.max_krylov_iters <<  ",  " << target_rank << ",  " << custom_rank << ",  " << residual_err_target <<  ",  " << residual_err_custom <<  ",  " << dur_rbki  << ",  " << dur_svd << ",\n";
end

function[Data_out] = data_preprocessing(Data_in, num_b_sizes, num_matmuls, numiters)
    
    Data_out = [];
    
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

function[] = plot_speed_iters(Data, num_b_sizes, num_matmuls, matmuls_to_display, display_mode)
    
    tiledlayout(1,2)
    legend_ctr = 1;
    marker_array = {'-o', '-diamond' '-s', '-^', '-v', '-+', '-*'};

    % Plot error vs #iters
    nexttile
    for i = 1 : num_b_sizes
        % number of iters performed = (target_rank / b_sz) * 2
        y = Data((i - 1) * num_matmuls + 1 : i * num_matmuls, 6); % or 6 for custom
        x = Data((i - 1) * num_matmuls + 1 : i * num_matmuls, 2);

        % Trim the data set
        x = x(1:matmuls_to_display, :);
        y = y(1:matmuls_to_display, :);

        if display_mode == 1
            % Remove even indices
            x(2:2:end, :) = [];
            y(2:2:end, :) = [];
        elseif display_mode == 2
            % Remove odd indices
            x(1:2:end, :) = [];
            y(1:2:end, :) = [];    
        end
        if i ~= 4 && i ~= 6
            plot(x, y, marker_array{i}, MarkerSize=12, LineWidth=2);
            set(gca,'YScale','log');
            hold on
        end
    end
    grid on
    xlabel('#GEMM(A)', 'FontWeight','bold') 
    ylabel('sqrt(||AV - SigmaU||^2 + ||A^TU-VSigma||^2)', 'FontWeight','bold')
    set(gca,'fontsize',14)

    % Plot error vs speedup over SVD
    nexttile
    xlim([0 inf])
    for i = 1 : num_b_sizes
        y = Data((i - 1) * num_matmuls + 1 : i * num_matmuls, 6); % or 6 for custom
        x = Data((i - 1) * num_matmuls + 1 : i * num_matmuls, 8) ./ Data((i - 1) * num_matmuls + 1 : i * num_matmuls, 7);
        
        % Trim the data set
        x = x(1:matmuls_to_display, :);
        y = y(1:matmuls_to_display, :);

        if display_mode == 1
            % Remove even indices
            x(2:2:end, :) = [];
            y(2:2:end, :) = [];
        elseif display_mode == 2
            % Remove odd indices
            x(1:2:end, :) = [];
            y(1:2:end, :) = [];    
        end

        if i ~= 4 && i ~= 6
            disp(i)
            plot(x, y, marker_array{i}, MarkerSize=12, LineWidth=2);
            set(gca,'YScale','log');
            set(gca,'XScale','log');
            legend_entries{legend_ctr} = ['b_{sz}=', num2str(Data(i * num_matmuls, 1))]; %#ok<AGROW>
            legend_ctr = legend_ctr + 1;
        end
        hold on
    end
    
    grid on
    xticks([0 50 100])
    xlabel('Speedup over SVD', 'FontWeight','bold') 
    lgd = legend(legend_entries, 'Location', 'southeast');
    fontsize(lgd, 14,'points')
    set(gca,'fontsize',14)
end