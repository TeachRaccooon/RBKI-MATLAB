function[U, S, VT] = RBKI(A, k, max_krylov_iters)

    [m, n] = size(A);
    iter_od = 0;
    iter_ev = 0;
    iter    = 0;
    Y_od = [];
    X_ev = [];
    
    Y_i = randn(n, k);
    X_i = A * Y_i;
    [X_ev, ~] = qr(X_i, 0);
    
    
    iter_od = iter_od + 1;

    for iter = 1:max_krylov_iters
        if (mod(iter, 2) ~= 0)
            Y_i = A' * X_i;
            if iter ~= 1
                R_i = Y_od' * Y_i;
                Y_i = Y_i - Y_od * R_i;
                Y_i = Y_i - Y_od *  Y_od' * Y_i;
            end
            [Y_i, R_ii] = qr(Y_i, 0);
            R = [R, R_i; zeros(size(R_ii, 1), size(R, 2)), R_ii]; %#ok<AGROW>
            Y_od = [Y_od, Y_i];                                   %#ok<AGROW>
            iter_ev = iter_ev + 1;
        else 
            X_i = A * Y_i;
            S_i = X_ev' * X_i;
            X_i = X_i - X_ev * S_i;
            X_i = X_i - X_ev * X_ev' * X_i;

            [X_i, S_ii] = qr(X_i, 0);
            S = [S, S_i; zeros(size(S_ii, 1), size(S, 2)), S_ii]; %#ok<AGROW>
            X_ev = [X_ev, X_i];                                   %#ok<AGROW>
            iter_od = iter_od + 1;
        end
    end

end