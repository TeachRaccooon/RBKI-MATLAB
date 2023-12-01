function[] = RBKI_benchmark()

end


function [] = call_RBKI(A, k, tol)
    A_cpy = A;
    [U1, Sigma1, V1] = RBKI_incremental_final(A, k, tol);
    [U2, Sigma2, V2] = svd(A_cpy);
    [U3, Sigma3, V3] = svdsketch(A_cpy, sqrt(tol));
    A_hat_rbki = U1(:, 1:k) * diag(Sigma1(1:k, 1)) * V1(:, 1:k)';
    A_hat_svd = U2(:, 1:k) * Sigma2(1:k, 1:k) * V2(:, 1:k)';
    fprintf("k: %d\n", k);
    fprintf("||S_svd  - S_rbki||_F/||S_svd||_F: %.20e\n",  norm(Sigma2(1:k, 1:k)  - diag(Sigma1(1:k, 1)), "fro") / norm(Sigma2(1:k, 1:k),  "fro"))
    fprintf("||S_svd  - S_rbki||_F/||S_svd||_F: %.20e\n",  norm(Sigma2(1:10, 1:10)  - diag(Sigma1(1:10, 1)), "fro") / norm(Sigma2(1:10, 1:10),  "fro"))
    %fprintf("||A      - A_rbki||_F/||A||_F: %.20e\n",      norm(A_cpy      - A_hat_rbki, "fro") / norm(A_cpy,      "fro"))
    fprintf("||A_svd  - A_rbki||_F/||A_svd||_F: %.20e\n",  norm(A_hat_svd  - A_hat_rbki, "fro") / norm(A_hat_svd,  "fro"))

    plot(Sigma1(1:k),            '-s', 'Color', 'red', "MarkerSize", 5,'LineWidth', 1)
    hold on
    plot(diag(Sigma2(1:k, 1:k)), '-o', 'Color', 'black', "MarkerSize", 5,'LineWidth', 1)
    grid on
    legend("RBKI", "SVD", "SVDsketch")
end

function[] =prepare_data()
    A = readmatrix("RBKI_test_matrrix_1.txt");
    % m = 957, n = 14079
    [m, n] = size(A);
    k = 1;


    while 1
    call_RBKI(A, k, tol);

    end
end

