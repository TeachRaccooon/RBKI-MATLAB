function[] = runtime_breakdown()
    
    % The first two entries in the dataset are: num_krylov_iters, b_sz
    Data_in = dlmread('../DATA_in/2024_10_re_running_all/RBKI_runtime_breakdown_m_10000_n_10000_k_start_16_k_stop_16_num_krylov_iters_start_16_num_krylov_iters_stop_16.txt');
    Data_out = [];

    [Data_in] = data_preprocessing_best(Data_in, 4, 2);

    for j = 1 : 4
        Data_out(j, 1) = 100 * Data_in(j, 3)                     /Data_in(j, 15); %#ok<AGROW> % Preallocation
        Data_out(j, 2) = 100 * Data_in(j, 4)                     /Data_in(j, 15); %#ok<AGROW> % SVD factors
        Data_out(j, 3) = 100 * Data_in(j, 5)                     /Data_in(j, 15); %#ok<AGROW> % UNGQR()
        Data_out(j, 4) = 100 * Data_in(j, 6)                     /Data_in(j, 15); %#ok<AGROW> % Reorth
        Data_out(j, 5) = 100 * Data_in(j, 7)                     /Data_in(j, 15); %#ok<AGROW> % QR
        Data_out(j, 6) = 100 * Data_in(j, 8)                     /Data_in(j, 15); %#ok<AGROW> % GEMM(A)
        %Data_out(j, 7) = 100 * Data_in(j, 9)                    /Data_in(j, 15); %#ok<AGROW> % EXCLUDE THIS - Main Loop  
        %Data_out(j, 7) = 100 * Data_in(j, 10)                    /Data_in(j, 15); %#ok<AGROW> % Sketching
        Data_out(j, 7) = 100 * (Data_in(j, 11) + Data_in(j, 12)) /Data_in(j, 15); %#ok<AGROW> % R_cpy + S_cpy
        %Data_out(j, 11) = 100 * Data_in(j, 13)                  /Data_in(j, 15); %#ok<AGROW> % Norm
        %Data_out(j, 5) = 100 * Data_in(j, 14)                   /Data_in(j, 15); %#ok<AGROW> % Rest
        
        Data_out(j, 8) = 100 * (Data_in(j, 10) + Data_in(j, 13) + Data_in(j, 14)) /Data_in(j, 15); %#ok<AGROW> % rest
        %Data_in(i, 15)
    end

    bplot = bar(Data_out,'stacked');
    bplot(1).FaceColor = 'black';
    bplot(2).FaceColor = 'yellow';
    bplot(3).FaceColor = 'green';
    bplot(4).FaceColor = 'magenta';
    bplot(5).FaceColor = '#EDB120';
    bplot(6).FaceColor = 'red';
    bplot(7).FaceColor = 'cyan';
    bplot(8).FaceColor = 'blue';
    
    bplot(1).FaceAlpha = 0.8;
    bplot(2).FaceAlpha = 0.8;
    bplot(3).FaceAlpha = 0.8;
    bplot(4).FaceAlpha = 0.8;
    bplot(5).FaceAlpha = 0.8;
    bplot(6).FaceAlpha = 0.8;
    bplot(7).FaceAlpha = 0.8;
    bplot(8).FaceAlpha = 0.8;
    
    lgd = legend('Data Alloc','SVD+Factors', 'UNGQR', 'Reorth', 'QR', 'GEMM(M)', 'Copy', 'Other');
    legend('Location','southeastoutside'); 
    set(gca,'XTickLabel',{'1', '4', '16', '64'});
    ylim([0 100]);
    ax = gca;
    ax.FontSize = 23; 
    lgd.FontSize = 20;

    title('ABRIK CPU', 'FontSize', 20);
    ylabel('Runtime %', 'FontSize', 20);
    xlabel('Block size', 'FontSize', 20);
end

function[Data_out] = data_preprocessing_best(Data_in, num_col_sizes, numiters)
    
    Data_out = [];
    i = 1;

    Data_out = [];
    while i < num_col_sizes * numiters
        best_speed = intmax;
        best_speed_idx = i;
        for j = 1:numiters
            if Data_in(i, 15) < best_speed
                best_speed = Data_in(i, 15);
                best_speed_idx = i;
            end
            i = i + 1;
        end
        Data_out = [Data_out; Data_in(best_speed_idx, :)]; %#ok<AGROW>
    end
end