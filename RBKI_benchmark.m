function[] = RBKI_benchmark()

    Data_out = prepare_data();
    writematrix(Data_out, 'DATA_in/RBKI_benchmark_out.txt')

end


function [] = call_RBKI(A, k, tol, numiters, Data_out)
    A_cpy = A;
    tic;
    [U1, Sigma1, V1] = RBKI_incremental_final(A, k, tol, numiters);
    t_rbki = toc;
    tic;
    [U2, Sigma2, V2] = svd(A_cpy);
    t_svd = toc;
    tic;
    [U3, Sigma3, V3] = RBKI_incremental_final(A, 1, tol, numiters);
    t_lanc = toc;

    
    err_rbki =  norm(Sigma2(1:k, 1:k)  - diag(Sigma1(1:k, 1)), "fro") / norm(Sigma2(1:k, 1:k),  "fro");
    err_lanc =  norm(Sigma2(1:1, 1:1)  - diag(Sigma3(1:1, 1)), "fro") / norm(Sigma2(1:1, 1:1),  "fro");
 
    fprintf("||S_svd  - S_rbki||_F/||S_svd||_F: %.20e\n", err_rbki);
    fprintf("||S_svd  - S_lanc||_F/||S_svd||_F: %.20e\n",  err_lanc);

    Data_out = [Data_out; k, numiters, err_rbki, err_lanc, t_rbki, t_svd, t_lanc];

end

function[Data_out] = prepare_data()
    A = readmatrix("DATA_in/RBKI_test_matrrix_1.txt");
    % m = 957, n = 14079
    [m, n] = size(A);
    tol = 2.5119e-14;

    b_sz = 2;
    b_sz_max = 4;
    numiters = 2;
    numiters_max = 64;
    numiters_start = numiters;

    Data_out = [];

    while b_sz <= b_sz_max
        while numiters < numiters_max
            call_RBKI(A, b_sz, tol, numiters, Data_out);
            numiters = numiters * 2;
        end
        numiters = numiters_start;
        b_sz = b_sz * 2;
    end
end

