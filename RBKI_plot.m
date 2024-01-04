%{
Below code plots the results of RBKI_speed_comp benchmark from RandLAPACK.
The following layout of an input data file is assumed by default:
Data column 1 (block size) varies from 2 to 'max_b_sz' in powers of 2.
Data column 2 (number of krylov iterations) varies from 2 to 
'max_krylov_iters' in powers of 2 per block size.
The file contains data for 5 input matrices from "Algorithm 971" paper, 
controlled by 'matrix_number' parameter.
%}

function[] = RBKI_plot()

    matrix_number = 1;
    max_b_sz = 2048;
    max_krylov_iters = 64;
    name = 'RBKI_benchmark_out.txt';
    Data = readmatrix(['DATA_in/' name]);
    
    plot_3d(Data, name, matrix_number, max_b_sz, max_krylov_iters);
    %plot_2d(Data, name, matrix_number, max_b_sz, max_krylov_iters);
end

function[] = plot_3d(Data, name, matrix_number, max_b_sz, max_krylov_iters)

    legend_entries = [];

    iter = 1;
    b_sz_power = log2(max_b_sz);
    krylov_iters_power = log2(max_krylov_iters);

    for i = (((matrix_number-1)*b_sz_power)+1):2:(((matrix_number-1)*b_sz_power)+b_sz_power)
        % RBKI singular values approximation accuracy
        y = Data((((i-1)*krylov_iters_power)+1):(i*krylov_iters_power), 3);
        % RBKI speedul over RSVD
        x = Data((((i-1)*krylov_iters_power)+1):(i * krylov_iters_power), krylov_iters_power) ./ Data((((i-1)*6)+1):(i * krylov_iters_power), 5);
        z = Data((((i-1)*krylov_iters_power)+1):(i * krylov_iters_power), 2);
    
        plot3(x, y, z, '-o', MarkerSize=10, LineWidth=2);
        set(gca,'YScale','log');
        hold on
        legend_entries{iter} = ['B_{sz}=', num2str(2^(i - ((matrix_number-1)*b_sz_power)))]; %#ok<AGROW>
        iter = iter + 1;
    end
    
    
    grid on

    title('RBKI benchmark')
    xlabel('Speedup over SVD') 
    ylabel('S.V. Accuracy')
    zlabel('#Krylov iters')
    legend(legend_entries);
    saveas(gcf, ['DATA_out/RBKI_benchmark_figures/Mat' int2str(matrix_number) '3d_' name '.fig'])
end

function[] = plot_2d(Data, name, matrix_number, max_b_sz, max_krylov_iters)

    legend_entries = [];

    iter = 1;
    b_sz_power = log2(max_b_sz)
    krylov_iters_power = log2(max_krylov_iters)

    for i = (((matrix_number-1)*b_sz_power)+1):2:(((matrix_number-1)*b_sz_power)+b_sz_power)
        % RBKI singular values approximation accuracy
        y = Data((((i-1)*krylov_iters_power)+1):(i*krylov_iters_power), 3);
        % RBKI speedul over RSVD
        x = Data((((i-1)*krylov_iters_power)+1):(i * krylov_iters_power), krylov_iters_power) ./ Data((((i-1)*krylov_iters_power)+1):(i * krylov_iters_power), 5);

        plot(x, y, '-o', MarkerSize=10, LineWidth=2);
        set(gca,'YScale','log');
        hold on
        legend_entries{iter} = ['B_{sz}=', num2str(2^(i - ((matrix_number-1)*b_sz_power)))]; %#ok<AGROW>
        iter = iter + 1;
    end
    
    
    grid on

    title('RBKI benchmark')
    xlabel('Speedup over SVD') 
    ylabel('S.V. Accuracy')
    legend(legend_entries);
    saveas(gcf, ['DATA_out/RBKI_benchmark_figures/large/Mat' int2str(matrix_number)  '_2d_' name '.fig'])
end