%Matrices from "Algorithm 971 paper"
function[] = gen_mat()

    m                 = 10000;
    n                 = 10000;
    rank              = 2000;
    low_rank          = 10;
    plotting_interval = 499;  %4999
    operation_mode = "generate";
    filepath = "Generated_matrices/test_mat_10k_rank_2k/";

    if operation_mode == "generate"
    
        U = randn(m, n);
        [U, ~] = qr(U, 0);
        V = randn(n, n);
        [V, ~] = qr(V, 0);
    
        Sigma = zeros(6, n);
        Sigma(1, :) = gen_mat_1(n);
        Sigma(2, :) = gen_mat_2(n, rank);
        Sigma(3, :) = gen_mat_3(n, rank);
        Sigma(4, :) = gen_mat_4(n, rank);
        Sigma(5, :) = gen_mat_5(n, rank);
        Sigma(6, :) = sort(abs(randn(1, n)),'descend');
   
        for i = 1:size(Sigma, 1)
            RBKI_test_mat_file = "RBKI_test_mat" + int2str(i) + ".txt";
            A_lowrank_mat_file = "A_lowrank_mat" + int2str(i) + ".txt";
            Spectrum_mat_file  = "Spectrum_mat"  + int2str(i) + ".txt";
        
            RBKI_test_mat_id = fopen(filepath + RBKI_test_mat_file, 'w');  
            A_lowrank_mat_id = fopen(filepath + A_lowrank_mat_file, 'w');  
            Spectrum_mat_id  = fopen(filepath + Spectrum_mat_file,  'w'); 
            fclose(RBKI_test_mat_id);
            fclose(A_lowrank_mat_id);
            fclose(Spectrum_mat_id);
        
            writematrix(U * diag(Sigma(i, :)) * V',                                        filepath + RBKI_test_mat_file,'delimiter',' ');
            writematrix(U(:, 1:low_rank) * diag(Sigma(i, 1:low_rank)) * V(:, 1:low_rank)', filepath + A_lowrank_mat_file,'delimiter',' ');
            writematrix(Sigma(i, :),                                                       filepath + Spectrum_mat_file ,'delimiter',' ');
            fprintf("Matrix " + int2str(i) + " processed\n");
        end
    elseif operation_mode == "plot"
        Sigma = zeros(6, n);
        for i = 1:6
            Sigma(i, :) = dlmread(filepath + "Spectrum_mat" + i + ".txt");
        end
    end
    plot(Sigma, n, plotting_interval, filepath) 
end

function [Sigma] = generator(n, k, type)
    switch type
        case 1
            Sigma = gen_mat_1(n);
        case 2
            Sigma = gen_mat_2(n, k);
        case 3
            Sigma = gen_mat_3(n, k);
        case 4
            Sigma = gen_mat_4(n, k);
        case 5
            Sigma = gen_mat_5(n, k);
        otherwise
            disp("Undefined type.")
    end
end

function[] = plot(Sigma, n, plotting_interval, filepath) 
    marker_array = {'-+', '-o', '-s', '-^', '-v', '-diamond', '-*'};
    x = 1:plotting_interval:n;
    semilogy(x, Sigma(1, 1:plotting_interval:end), marker_array{1}, MarkerSize=18, LineWidth=1.8)
    hold on
    semilogy(x, Sigma(2, 1:plotting_interval:end), marker_array{2}, MarkerSize=18, LineWidth=1.8)
    hold on
    semilogy(x, Sigma(3, 1:plotting_interval:end), marker_array{3}, MarkerSize=28, LineWidth=1.8)
    hold on
    semilogy(x, Sigma(4, 1:plotting_interval:end), marker_array{4}, MarkerSize=18, LineWidth=1.8)
    hold on
    semilogy(x, Sigma(5, 1:plotting_interval:end), marker_array{5}, MarkerSize=18, LineWidth=1.8)
    hold on
    semilogy(x, Sigma(6, 1:plotting_interval:end), marker_array{6}, MarkerSize=18, LineWidth=1.8)
    
    grid on
    xlabel('k', 'FontSize', 20);
    ylabel('sigma[k]', 'FontSize', 20);
    lgd = legend('Mat 1', 'Mat 2', 'Mat 3', 'Mat 4', 'Mat 5', 'Mat 6');
    lgd.FontSize = 20;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;

    saveas(gcf, filepath + 'specta_plot.fig')
end

function[Sigma] = gen_mat_1(n)
    Sigma = zeros(1, n);
    for j = 1:n
        Sigma(1, j) = 1/j;
    end
end

function[Sigma] = gen_mat_2(n, k)
    Sigma = zeros(1, n);
    Sigma(1, 1) = 1;
    for j = 2:k
        Sigma(1, j) = 2 * 10^(-5);
    end
    for j = (k+1):n
        Sigma(1, j) = 10^(-5) * (k + 1) / j;
    end
end

function[Sigma] = gen_mat_3(n, k)
    Sigma = zeros(1, n);
    for j = 1:k
        Sigma(1, j) = 10^(-5 * (j - 1) / (k - 1));
    end
    for j = (k+1):n
        Sigma(1, j) = 10^(-5) * (k + 1) / j;
    end
end

function[Sigma] = gen_mat_4(n, k)
    Sigma = zeros(1, n);
    for j = 1:k
        Sigma(1, j) = 10^(-5 * (j - 1) / (k - 1));
    end
    Sigma(1, j) = 10^(-5);
    for j = (k+2):n
        Sigma(1, j) = 0;
    end
end

function[Sigma] = gen_mat_5(n, k)
    Sigma = zeros(1, n);
    for j = 1:k
        Sigma(1, j) = 10^(-5) + (1 - 10^(-5)) * (k - j) / (k - 1);
    end
    for j = (k+1):n
        Sigma(1, j) = 10^(-5) * sqrt((k + 1) / j);
    end
end


