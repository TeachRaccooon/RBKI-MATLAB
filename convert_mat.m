% Fortran is pain, so it is easier to just remake my data scripts
% This functions transforms the generated matrix into a one-dimensional
% column-major format, easily dijestible by fortran.
function[] = convert_mat()
    for i = 1:5
        Data = readmatrix(['DATA_out/test_mat_small/RBKI_test_mat' num2str(i) '.txt']);
        [m, n] = size(Data);

        writematrix(reshape(Data, [m*n, 1]), ['DATA_out/fortran_test_mat_small/RBKI_test_mat' num2str(i) '.txt'],'delimiter',' ')
    end
end