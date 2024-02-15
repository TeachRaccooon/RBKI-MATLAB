%Matrices from "Algorithm 971 paper"
function[] = gem_mat()

    m = 10000;
    n = 10000;
    k = 2000;

    U = randn(m, n);
    [U, ~] = qr(U, 0);
    V = randn(n, n);
    [V, ~] = qr(V, 0);


    [Sigma1] = gen_mat_1(n);
    [Sigma2] = gen_mat_2(n, k);
    [Sigma3] = gen_mat_3(n, k);
    [Sigma4] = gen_mat_4(n, k);
    [Sigma5] = gen_mat_5(n, k);

    writematrix(U * diag(Sigma1) * V', 'DATA_in/test_mat_10k_rank_2k/RBKI_test_mat1.txt','delimiter',' ');
    writematrix(Sigma1               , 'DATA_in/test_mat_10k_rank_2k/spectrum_mat1.txt','delimiter',' ');
    fprintf("Done with Mat1\n");
    writematrix(U * diag(Sigma2) * V', 'DATA_in/test_mat_10k_rank_2k/RBKI_test_mat2.txt','delimiter',' ');
    writematrix(Sigma2               , 'DATA_in/test_mat_10k_rank_2k/spectrum_mat2.txt','delimiter',' ');
    fprintf("Done with Mat2\n");
    writematrix(U * diag(Sigma3) * V', 'DATA_in/test_mat_10k_rank_2k/RBKI_test_mat3.txt','delimiter',' ');
    writematrix(Sigma3               , 'DATA_in/test_mat_10k_rank_2k/spectrum_mat3.txt','delimiter',' ');
    fprintf("Done with Mat3\n");
    writematrix(U * diag(Sigma4) * V', 'DATA_in/test_mat_10k_rank_2k/RBKI_test_mat4.txt','delimiter',' ');
    writematrix(Sigma4               , 'DATA_in/test_mat_10k_rank_2k/spectrum_mat4.txt','delimiter',' ');
    fprintf("Done with Mat4\n");
    writematrix(U * diag(Sigma5) * V', 'DATA_in/test_mat_10k_rank_2k/RBKI_test_mat5.txt','delimiter',' ');
    writematrix(Sigma5               , 'DATA_in/test_mat_10k_rank_2k/spectrum_mat5.txt','delimiter',' ');
    fprintf("Done with Mat5\n");

    x = 1:n;
    semilogy(x, Sigma1, '-o', MarkerSize=3, LineWidth=2)
    hold on
    semilogy(x, Sigma2, '-o', MarkerSize=3, LineWidth=2)
    hold on
    semilogy(x, Sigma3, '-o', MarkerSize=3, LineWidth=2)
    hold on
    semilogy(x, Sigma4, '-o', MarkerSize=3, LineWidth=2)
    hold on
    semilogy(x, Sigma5, '-o', MarkerSize=3, LineWidth=2)
    
    
    legend('1', '2', '3', '4', '5')
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


