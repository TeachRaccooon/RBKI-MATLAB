function[] = runtime_breakdown()
    
    % The first two entries in the dataset are: num_krylov_iters, b_sz
    Data_in = dlmread('DATA_in/2024_04_08_edelman_mat/Edelman_small_runtime_breakdown.txt');
    Data_out = [];

    lines_we_want = [3, 5, 7];

    Data_in(7, 1)

    for j = 1 : 3

        Data_out(j, 1) = 100 * Data_in(lines_we_want(1, j), 3)                  /Data_in(lines_we_want(1, j), 15); %#ok<AGROW> % Preallocation
        Data_out(j, 2) = 100 * Data_in(lines_we_want(1, j), 4)                  /Data_in(lines_we_want(1, j), 15); %#ok<AGROW> % SVD factors
        Data_out(j, 3) = 100 * Data_in(lines_we_want(1, j), 5)                  /Data_in(lines_we_want(1, j), 15); %#ok<AGROW> % UNGQR()
        Data_out(j, 4) = 100 * Data_in(lines_we_want(1, j), 6)                  /Data_in(lines_we_want(1, j), 15); %#ok<AGROW> % Reorth
        Data_out(j, 5) = 100 * Data_in(lines_we_want(1, j), 7)                  /Data_in(lines_we_want(1, j), 15); %#ok<AGROW> % QR
        Data_out(j, 6) = 100 * Data_in(lines_we_want(1, j), 8)                  /Data_in(lines_we_want(1, j), 15); %#ok<AGROW> % GEMM(A)
        %Data_out(j, 7) = 100 * Data_in(lines_we_want(1, j), 9)                  /Data_in(lines_we_want(1, j), 15); %#ok<AGROW> % EXCLUDE THIS - Main Loop  
        Data_out(j, 7) = 100 * Data_in(lines_we_want(1, j), 10)                  /Data_in(lines_we_want(1, j), 15); %#ok<AGROW> % Sketching
        %Data_out(j, 9) = 100 * Data_in(lines_we_want(1, j), 11)                  /Data_in(lines_we_want(1, j), 15); %#ok<AGROW> % R_cpy
        %Data_out(j, 10) = 100 * Data_in(lines_we_want(1, j), 12)                /Data_in(lines_we_want(1, j), 15); %#ok<AGROW> % S_cpy
        %Data_out(j, 11) = 100 * Data_in(lines_we_want(1, j), 13)                /Data_in(lines_we_want(1, j), 15); %#ok<AGROW> % Norm
        %Data_out(j, 5) = 100 * Data_in(lines_we_want(1, j), 14)                /Data_in(lines_we_want(1, j), 15); %#ok<AGROW> % Rest
        
        Data_out(j, 8) = 100 * (Data_in(lines_we_want(1, j), 11) + Data_in(lines_we_want(1, j), 12) + Data_in(lines_we_want(1, j), 13) + Data_in(lines_we_want(1, j), 14)) /Data_in(lines_we_want(1, j), 15); %#ok<AGROW> % rest
        %Data_in(i, 15)
    end

    Data_out

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
    
    lgd = legend('Data Alloc','SVD+Factors', 'UNGQR', 'Reorth', 'QR', 'GEMM(A)', 'Sketching', 'Other')
    legend('Location','northeastoutside'); 
    set(gca,'XTickLabel',{50', '75', '120'});
    ylim([0 100]);
    ax = gca;
    ax.FontSize = 23; 
    lgd.FontSize = 15;

    %title('Mat5 b_sz 128 runtime breakdown')
    %saveas(gcf,'DATA_out/test_mat_100k_rank_20k/Mat1_b_sz_128_breakdown.jpg')
end