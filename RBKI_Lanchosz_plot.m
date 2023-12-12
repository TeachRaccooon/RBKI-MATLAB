function[] = RBKI_Lanchosz_plot()

    name = 'RBKI_speed_comp_m_5000_n_5000_k_start_2_k_stop_256_num_krylov_iters_start_2_num_krylov_iters_stop_64.dat';
    Data = readmatrix(['DATA_in/' name]);
    
    %plot_3d(Data, name);
    plot_2d(Data, name);
end

function[] = plot_3d(Data, name)

    legend_entries = [];

    iter = 1;
    for i = 1:2:9
        % RBKI_Lan singular values approximation accuracy
        y = Data((((i-1)*6)+1):(i*6), 4);
        % RBKI_Lan speedul over RSVD
        x = Data((((i-1)*6)+1):(i * 6), 6) ./ Data((((i-1)*6)+1):(i * 6), 7);
        z = Data((((i-1)*6)+1):(i * 6), 2);
    
        plot3(x, y, z, '-o', MarkerSize=10, LineWidth=2);
        set(gca,'YScale','log');
        hold on
        legend_entries{iter} = ['B_{sz}=', num2str(2^i)]; %#ok<AGROW>
        iter = iter + 1;
    end
    
    
    grid on

    title('RBKI Lanchosz benchmark')
    xlabel('Speedup over SVD') 
    ylabel('S.V. Accuracy')
    zlabel('#Krylov iters')
    legend(legend_entries);
    saveas(gcf, ['DATA_out/RBKI_benchmark_figures/3d_Lan_' name '.fig'])
end

function[] = plot_2d(Data, name)

    legend_entries = [];

    iter = 1;
    for i = 1:2:9
        % RBKI_Lan singular values approximation accuracy
        y = Data((((i-1)*6)+1):(i*6), 4);
        % RBKI_Lan speedul over RSVD
        x = Data((((i-1)*6)+1):(i * 6), 6) ./ Data((((i-1)*6)+1):(i * 6), 7);

        plot(x, y, '-o', MarkerSize=10, LineWidth=2);
        set(gca,'YScale','log');
        xlim([0 250]);
        hold on
        legend_entries{iter} = ['B_{sz}=', num2str(2^i)]; %#ok<AGROW>
        iter = iter + 1;
    end
    
    
    grid on

    title('RBKI Lanchosz benchmark')
    xlabel('Speedup over SVD') 
    ylabel('S.V. Accuracy')
    legend(legend_entries);
    saveas(gcf, ['DATA_out/RBKI_benchmark_figures/2d_Lan_' name '.fig'])
end