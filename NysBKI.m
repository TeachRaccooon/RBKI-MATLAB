function[U, Lambda] = NysBKI(A, k, tol)
    [~, n] = size(A);
    norm_A = norm(A, 'fro');
    sq_tol = tol^2;
    [X_i, ~] = qr(randn(n, k), 0);
    Y_i = A * X_i;
    S_i = [];
    Y_full = Y_i;
    X_full = X_i;
    iter = 1;
    while 1
        iter
        %X_prev = X_full(:, (k*(iter-1)+1):end);
        X_i = Y_i + tol * X_i;
        if iter == 1
            R_i = X_full' * X_i;
        else
            R_i = [zeros(size(S_i, 1) - 2 * k), X_full(:, (k*(iter-2)+1):(k*(iter-1))), X_full(:, (k*(iter-1)+1):end)]' * X_i;
        end
        X_i = X_i - X_full * X_full' * X_i;
        X_i = X_i - X_full * X_full' * X_i;
        [X_i, R_ii] = qr(X_i, 0);
        X_full = [X_full X_i]; %#ok<AGROW>

        % Termination criteria
        try chol([S_i R_i]);
            C = chol([S_i R_i]);
        catch
            break;
        end

        S_i = [S_i R_i; zeros(size(R_ii, 1), size(S_i, 2)) R_ii]; %#ok<AGROW>
        Z = S_i * inv(C);
        [U_hat, Sigma, ~] = svd(Z, 0);
        U = X_full * U_hat;
        for i = 1:size(Sigma, 1)
            Lambda(i, i) = max(Sigma(i, i)^2 - eps, 0); %#ok<AGROW>
        end
        %Lambda = diag(Sigma)^2 - tol * eye(size(Sigma, 1));
        Y_i = A * X_i;
        Y_full = [Y_full Y_i];  %#ok<AGROW>
        E = Y_full * U_hat - U * Sigma;
        norm(E, 'fro')
        iter = iter + 1;
    end
end