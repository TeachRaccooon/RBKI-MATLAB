% This code is temporary bc idk how to properly incorporate the fact that
% some of the data for block sizes is missing
function[] = RBKI_benchmark_Riley_analysis()
    Data_in = readmatrix("DATA_in/2024_04_08_edelman_mat/EDELMAN_MAT_speed_comp_m_41505_n_81700_b_sz_start_1_b_sz_stop_512_num_matmuls_start_2_num_matmuls_stop_75.txt");
    
    for i = 1:size(Data_in, 1)
        % Mat 1 SVD:        1860675838
        % Mat 2 SVD:        1966542609
        % Mat 3 SVD:        1824193427
        % Mat 4 SVD:        2629361560
        % Mat 5 SVD:        2587508196
        % Edelman Matrix 1: 517627648
        % Edelman Matrix 2: 2576854142
        Data_in(i, 8) = 517627648;
    end
    size(Data_in)
    
    min_b_sz = 1;
    max_b_sz = 512;

    max_matmuls = 75;
    min_matmuls = 2;

    % We don't have b_sz 4 in Mat 1-5 dataset, so need to add -1 at the end of the line below.
    num_b_sizes = (log2(max_b_sz) - log2(min_b_sz) + 1)
    num_matmuls = max_matmuls - min_matmuls + 1


   % for i = 2:1436
   %     if Data_in(i, 1) ~= Data_in(i - 1, 1)
   %     fprintf("%d\n", i)
   %     end
   % end


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
    [m, n] = size(Data);
    tiledlayout(1,2)
    legend_ctr = 1;
    marker_array = {'-o', '-diamond' '-s', '-^', '-v', '-+', '-*', '->', '-pentagram' '<-'};

    % Plot error vs #iters
    nexttile
    for i = 1 : num_b_sizes
        % number of iters performed = (target_rank / b_sz) * 2
        y = Data((i - 1) * num_matmuls + 1 : i * num_matmuls, 6) / sqrt(m + n); % or 6 for custom
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
    for i = 1 : num_b_sizes
        y = Data((i - 1) * num_matmuls + 1 : i * num_matmuls, 6) / sqrt(m + n); % or 6 for custom
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

        %if i ~= 4 && i ~= 6
            disp(i)
            plot(x, y, marker_array{i}, MarkerSize=12, LineWidth=2);
            set(gca,'YScale','log');
            set(gca,'XScale','log');
            legend_entries{legend_ctr} = ['b_{sz}=', num2str(Data(i * num_matmuls, 1))]; %#ok<AGROW>
            legend_ctr = legend_ctr + 1;
        %end
        hold on
    end
    
    grid on
    xlabel('Speedup over SVD', 'FontWeight','bold') 
    lgd = legend(legend_entries, 'Location', 'southeast');
    fontsize(lgd, 14,'points')
    set(gca,'fontsize',14)
end

