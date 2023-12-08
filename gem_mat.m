%Matrices from "Algorithm 971 paper"
function[] = gem_mat()

    m = 5000;
    n = 5000;
    k = 100;

    [A1, Sigma1] = gen_mat_1(m, n);
    [A2, Sigma2] = gen_mat_2(m, n, k);
    [A3, Sigma3] = gen_mat_3(m, n, k);
    [A4, Sigma4] = gen_mat_4(m, n, k);
    [A5, Sigma5] = gen_mat_5(m, n, k);

    writematrix(A1, 'DATA_out/test_mat_med/RBKI_test_mat1.txt','delimiter',' ')
    writematrix(A2, 'DATA_out/test_mat_med/RBKI_test_mat2.txt','delimiter',' ');
    writematrix(A3, 'DATA_out/test_mat_med/RBKI_test_mat3.txt','delimiter',' ');
    writematrix(A4, 'DATA_out/test_mat_med/RBKI_test_mat4.txt','delimiter',' ');
    writematrix(A5, 'DATA_out/test_mat_med/RBKI_test_mat5.txt','delimiter',' ');
    
    x = 1:n;
    semilogy(x, diag(Sigma1), '-o', MarkerSize=3, LineWidth=2)
    hold on
    semilogy(x, diag(Sigma2), '-o', MarkerSize=3, LineWidth=2)
    hold on
    semilogy(x, diag(Sigma3), '-o', MarkerSize=3, LineWidth=2)
    hold on
    semilogy(x, diag(Sigma4), '-o', MarkerSize=3, LineWidth=2)
    hold on
    semilogy(x, diag(Sigma5), '-o', MarkerSize=3, LineWidth=2)

    legend('1', '2', '3', '4', '5')
end


function[A, Sigma] = gen_mat_1(m, n)
    [U, ~] = qr(randn(m, n), 0);
    [V, ~] = qr(randn(n, n), 0);
    Sigma = zeros(n, n);
    for j = 1:n
        Sigma(j, j) = 1/j;
    end
    A = U * Sigma * V';
end

function[A, Sigma] = gen_mat_2(m, n, k)
    [U, ~] = qr(randn(m, n), 0);
    [V, ~] = qr(randn(n, n), 0);
    Sigma = zeros(n, n);
    Sigma(1, 1) = 1;
    for j = 2:k
        Sigma(j, j) = 2 * 10^(-5);
    end
    for j = (k+1):n
        Sigma(j, j) = 10^(-5) * (k + 1) / j;
    end
    A = U * Sigma * V';
end

function[A, Sigma] = gen_mat_3(m, n, k)
    [U, ~] = qr(randn(m, n), 0);
    [V, ~] = qr(randn(n, n), 0);
    Sigma = zeros(n, n);
    for j = 1:k
        Sigma(j, j) = 10^(-5 * (j - 1) / (k - 1));
    end
    for j = (k+1):n
        Sigma(j, j) = 10^(-5) * (k + 1) / j;
    end
    A = U * Sigma * V';
end

function[A, Sigma] = gen_mat_4(m, n, k)
    [U, ~] = qr(randn(m, n), 0);
    [V, ~] = qr(randn(n, n), 0);
    Sigma = zeros(n, n);
    for j = 1:k
        Sigma(j, j) = 10^(-5 * (j - 1) / (k - 1));
    end
    Sigma(j, j) = 10^(-5);
    for j = (k+2):n
        Sigma(j, j) = 0;
    end
    A = U * Sigma * V';
end

function[A, Sigma] = gen_mat_5(m, n, k)
    [U, ~] = qr(randn(m, n), 0);
    [V, ~] = qr(randn(n, n), 0);
    Sigma = zeros(n, n);
    for j = 1:k
        Sigma(j, j) = 10^(-5) + (1 - 10^(-5)) * (k - j) / (k - 1);
    end
    for j = (k+1):n
        Sigma(j, j) = 10^(-5) * sqrt((k + 1) / j);
    end
    A = U * Sigma * V';
end


