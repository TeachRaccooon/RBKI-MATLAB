function[] = convert_mat_fortran()

    Data_in = readmatrix("RBKI_test_mat1");
    [m, n] = size(Data_in);
    
    Data_out = [];
    for i = 1:n
        Data_out = [Data_out; Data_in(:, i)];
    end
    writematrix(Data_out, 'DATA_out/RBKI_test_mat1.txt','delimiter',' ');

end