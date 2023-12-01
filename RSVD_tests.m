
function[] = RSVD_tests()
    NysBKI_test()
    %RBKI_DATASET_test()
    %RBKI_random_test()
    %RBKI_hard_test
end
function [] = call_RBKI(A, k, tol)
    A_cpy = A;
    [U1, Sigma1, V1] = RBKI_incremental_final(A, k, tol);
    [U2, Sigma2, V2] = svd(A_cpy);
    [U3, Sigma3, V3] = svdsketch(A_cpy, sqrt(tol));
    A_hat_rbki = U1(:, 1:k) * diag(Sigma1(1:k, 1)) * V1(:, 1:k)';
    A_hat_svd = U2(:, 1:k) * Sigma2(1:k, 1:k) * V2(:, 1:k)';
    A_hat_rsvd = U3(:, 1:k) * Sigma3(1:k, 1:k) * V3(:, 1:k)';
    fprintf("k: %d\n", k);
    fprintf("||S_svd  - S_rbki||_F/||S_svd||_F: %.20e\n",  norm(Sigma2(1:k, 1:k)  - diag(Sigma1(1:k, 1)), "fro") / norm(Sigma2(1:k, 1:k),  "fro"))
    fprintf("||S_svd  - S_rbki||_F/||S_svd||_F: %.20e\n",  norm(Sigma2(1:10, 1:10)  - diag(Sigma1(1:10, 1)), "fro") / norm(Sigma2(1:10, 1:10),  "fro"))
    %fprintf("||A      - A_rbki||_F/||A||_F: %.20e\n",      norm(A_cpy      - A_hat_rbki, "fro") / norm(A_cpy,      "fro"))
    fprintf("||A_svd  - A_rbki||_F/||A_svd||_F: %.20e\n",  norm(A_hat_svd  - A_hat_rbki, "fro") / norm(A_hat_svd,  "fro"))
    fprintf("||A_rsvd - A_rbki||_F/||A_rsvd||_F: %.20e\n", norm(A_hat_rsvd - A_hat_rbki, "fro") / norm(A_hat_rsvd, "fro"))
    plot(Sigma1(1:k),            '-s', 'Color', 'red', "MarkerSize", 5,'LineWidth', 1)
    hold on
    plot(diag(Sigma2(1:k, 1:k)), '-o', 'Color', 'black', "MarkerSize", 5,'LineWidth', 1)
    hold on
    plot(diag(Sigma3(1:k, 1:k)), '-^', 'Color', 'blue', "MarkerSize", 5,'LineWidth', 1)
    grid on
    legend("RBKI", "SVD", "SVDsketch")
end
function [] = call_NysBKI(A, k, tol)
    A_cpy = A;
    [U1, Lambda1] = NysBKI(A, k, tol);
    [U2, Lambda2] = eig(A_cpy);
    A_hat_rbki = U1(:, 1:k) * Lambda1(1:k, 1:k) * U1(:, 1:k)';
    A_hat_evd  = U2(:, 1:k) * Lambda2(1:k, 1:k) * U2(:, 1:k)';
    fprintf("k: %d\n", k);
    fprintf("||L_evd  - L_nysbki||_F/||L_evd||_F: %.20e\n",  norm(Lambda2(1:k, 1:k)  - Lambda1(1:k, 1:k), "fro") / norm(Lambda2(1:k, 1:k),  "fro"))
    fprintf("||A_evd  - A_nysbki||_F/||A_evd||_F: %.20e\n",  norm(A_hat_evd  - A_hat_rbki, "fro") / norm(A_hat_evd,  "fro"))
    plot(diag(Lambda1(1:k, 1:k)), '-s', 'Color', 'red', "MarkerSize", 5,'LineWidth', 1)
    hold on
    plot(diag(Lambda2(1:k, 1:k)), '-o', 'Color', 'black', "MarkerSize", 5,'LineWidth', 1)
    grid on
    legend("NysBKI", "EVD")
end
function[] = RBKI_DATASET_test()
    fprintf("/--------------------------------------------------------/\n")
    A = readmatrix("RBKI_test_matrrix_1.txt");
    %A = A';
    % m = 957, n = 14079
    [m, n] = size(A);
    A = randn(m, n);
    tol = 2.5119e-14;
    k = 1;
    call_RBKI(A, k, tol);
    fprintf("/--------------------------------------------------------/\n")
end
function[] = RBKI_random_test()
    fprintf("/--------------------------------------------------------/\n")
    m = 2000;
    n = 2000;
    tol = 2.5119e-14;
    k = 999;
    A = randn(m, n);
    A = diag(diag(A));
    call_RBKI(A, k, tol);
    fprintf("/--------------------------------------------------------/\n")
end
function[] = RBKI_hard_test()
    fprintf("/--------------------------------------------------------/\n")
    m = 2000;
    n = 2000;
    tol = 2.5119e-14;
    k = 35;
    A = zeros(1, n);
    A(:, 1:9) = [10, 9, 8, 7, 6, 5, 4, 3, 2];
    for i = 10:k
        A(1, i) = A(1, i-1) - (1 / n);
    end
    A = diag(A);
    call_RBKI(A, k, tol);
    fprintf("/--------------------------------------------------------/\n")
end
function[] = NysBKI_test()
    fprintf("/--------------------------------------------------------/\n")
    m = 5;
    tol = 2.5119e-14;
    k = 5;
    A = zeros(1, m);
    A(:, 1:m) = [10, 9, 8, 7, 6];
    A = diag(A(1, 1:5));
    [B, ~] = qr(randn(m, m), 0);
    A = B * B' * A;

    call_NysBKI(A, k, tol);
    fprintf("/--------------------------------------------------------/\n")
end