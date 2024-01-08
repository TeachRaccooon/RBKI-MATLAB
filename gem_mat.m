%Matrices from "Algorithm 971 paper"
function[] = gem_mat()

    m = 100000;
    n = 100000; 
    k = 2000;

    U = randn(m, n);
    U = orth(U, 0);
    writematrix(U, 'DATA_in/test_mat_100k/U.txt','delimiter',' ');
    clear U;

    V = randn(n, n);
    V = orth(V, 0);
    writematrix(V', 'DATA_in/test_mat_100k/V.txt','delimiter',' ');
    clear V;
    
    for type = 1:5
        [Sigma] = generator(n, k, type);
        writematrix(Sigma, 'DATA_in/test_mat_100k/Sigma.txt','delimiter',' ');
        clear Sigma;
        % Goal is to assemble 100k by 100k matrix in 10k by 100k pieces. 
        % Load 10k by 50k piece of U
        % Load 50k values of Sigma
        % load 50 k by 50 k piece of V
    
        mat_sz        = 10;
        % Number of wors in a single stored block
        % Note that .* works only with 1-d vectors
        block_sz      = 1; 
        
        for start_row_idx = 1:100000
            %start_row_idx
            % Fill b_sz by mat_sz / 2 block with data;
            %U_small = readmatrix('DATA_in/test_mat_100k/U.txt', 'Range',     [start_row_idx, 1, start_row_idx + block_sz - 1, mat_sz]);
            %V_small = readmatrix('DATA_in/test_mat_100k/V.txt', 'Range',     [1, 1, mat_sz, mat_sz / 2]);
            %S_small = readmatrix('DATA_in/test_mat_100k/Sigma.txt', 'Range', [1, 1, mat_sz, mat_sz]);
            %A = U_small .* S_small * V_small;
            A = readmatrix('DATA_in/test_mat_100k/U.txt', 'Range',     [start_row_idx, 1, start_row_idx + block_sz - 1, mat_sz]) .* readmatrix('DATA_in/test_mat_100k/Sigma.txt', 'Range', [1, 1, mat_sz, mat_sz]) * readmatrix('DATA_in/test_mat_100k/V.txt', 'Range',     [1, 1, mat_sz, mat_sz / 2]);
        
            % Fill b_sz by mat_sz block with data;
            %U_small = readmatrix('DATA_in/test_mat_100k/U.txt', 'Range',     [start_row_idx, 1, start_row_idx + block_sz - 1, mat_sz]);
            %V_small = readmatrix('DATA_in/test_mat_100k/V.txt', 'Range',     [1, (mat_sz / 2) + 1, mat_sz, mat_sz]);
            %S_small = readmatrix('DATA_in/test_mat_100k/Sigma.txt', 'Range', [1, block_sz, mat_sz, mat_sz]);
            %A = U_small .* S_small * V_small;
            A = [A, readmatrix('DATA_in/test_mat_100k/U.txt', 'Range',     [start_row_idx, 1, start_row_idx + block_sz - 1, mat_sz]) .* readmatrix('DATA_in/test_mat_100k/Sigma.txt', 'Range', [1, block_sz, mat_sz, mat_sz]) * readmatrix('DATA_in/test_mat_100k/V.txt', 'Range',     [1, (mat_sz / 2) + 1, mat_sz, mat_sz])];
            
            writematrix(A, ['DATA_in/test_mat_100k/RBKI_test_mat' int2str(type) '.txt'], 'delimiter', ' ', 'WriteMode','append');
            clear A;
        end
    end

    %[Sigma2] = gen_mat_2(n, k);
    %[Sigma3] = gen_mat_3(n, k);
    %[Sigma4] = gen_mat_4(n, k);
    %[Sigma5] = gen_mat_5(n, k);

    %writematrix(A1, 'DATA_out/test_mat_large/RBKI_test_mat1.txt','delimiter',' ');
    %{
    writematrix(A2, 'DATA_out/test_mat_large/RBKI_test_mat2.txt','delimiter',' ');
    writematrix(A3, 'DATA_out/test_mat_large/RBKI_test_mat3.txt','delimiter',' ');
    writematrix(A4, 'DATA_out/test_mat_large/RBKI_test_mat4.txt','delimiter',' ');
    writematrix(A5, 'DATA_out/test_mat_large/RBKI_test_mat5.txt','delimiter',' ');
    
    x = 1:n;
    semilogy(x, diag(Sigma1), '-o', MarkerSize=3, LineWidth=2)
    hold on
    semilogy(x, diag(Sigma2), '-o', MarkerSize=10, LineWidth=2)
    hold on
    semilogy(x, diag(Sigma3), '-o', MarkerSize=3, LineWidth=2)
    hold on
    semilogy(x, diag(Sigma4), '-o', MarkerSize=3, LineWidth=2)
    hold on
    semilogy(x, diag(Sigma5), '-o', MarkerSize=3, LineWidth=2)
    
    
    legend('1', '2', '3', '4', '5')
    %}
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


