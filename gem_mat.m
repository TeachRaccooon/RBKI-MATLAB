%Matrices from "Algorithm 971 paper"
function[] = gem_mat()

    m = 5;
    n = 5;
    k = 2;

    U = randn(m, n);
    [U, ~] = qr(U, 0);
    writematrix(U, 'DATA_in/test_mat_100k/U.txt','delimiter',' ');
    U1 = U;
    clear U;

    V = randn(n, n);
    [V, ~] = qr(V, 0);
    writematrix(V, 'DATA_in/test_mat_100k/V.txt','delimiter',' ');
    V1 = V;
    clear V;
    
    [Sigma] = gen_mat_1(n);
    writematrix(Sigma, 'DATA_in/test_mat_100k/Sigma.txt','delimiter',' ');
    S1 = Sigma;
    clear Sigma1;

    A = U1 .* S1 * V1

    U_small = readmatrix('DATA_in/test_mat_100k/U.txt', 'Range',     [1, 1, 2, 10])
    V_small = readmatrix('DATA_in/test_mat_100k/V.txt', 'Range',     [1, 1, 10, 2])
    S_small = readmatrix('DATA_in/test_mat_100k/Sigma.txt', 'Range', [1, 1, 1, 2])


    writematrix(A_small, 'DATA_in/test_mat_100k/A_small.txt','delimiter',' ');



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


