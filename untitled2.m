

A = randn(1000, 1000);
B = A;
tic;
svd(A);
t_svd = toc;

tic;
svds(B, 20);
t_svds = toc;

fprintf("Time SVD: %f\n", t_svd);
fprintf("Time SVDS: %f\n", t_svds);


