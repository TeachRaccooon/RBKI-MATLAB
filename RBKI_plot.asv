function[] = RBKI_plot()

    Data = readmatrix('DATA_in/RBKI_try.txt');
    legend_entries = [];

    iter = 1;
    for i = 1:2:9
        % RBKI singular values approximation accuracy
        y = Data((((i-1)*6)+1):(i*6), 3);
        % RBKI speedul over RSVD
        x = Data((((i-1)*6)+1):(i * 6), 6) ./ Data((((i-1)*6)+1):(i * 6), 5);
        z = Data((((i-1)*6)+1):(i * 6), 2);
    
        plot3(x, y, z, '-o', MarkerSize=10, LineWidth=2);
        set(gca,'YScale','log');
        hold on
        legend_entries{iter} = ['B_{sz}=', num2str(2^i)]; %#ok<AGROW>
        iter = iter + 1;
    end
    
    
    grid on

    xlabel('Speedup over SVD') 
    ylabel('S.V. Accuracy')
    zlabel('#Krylov iters')
    legend(legend_entries);
end