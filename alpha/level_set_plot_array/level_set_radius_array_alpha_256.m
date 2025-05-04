function [] = level_set_radius_array_alpha(R0, m, alpha, Nx)
    indir="/scratch/xpz5km/Cahn_Hilliard_Model/julia_out/critical_radius_alpha"
    folder = "plots";
    total_time=10;
    everyR=10;

    epsilon = m * (1 / Nx) / (2 * sqrt(2) * atanh(0.9))
    epsilon_name = sprintf('%.5g', epsilon)
    
    R0_name = sprintf('%.5g', R0)
    alpha_name = sprintf('%.1f', alpha)


    [rr,tt] = level_set_plot_alpha_256(2.5e-5, indir, total_time, everyR, epsilon_name, R0_name, folder, alpha_name, Nx);
    R0_vector = repmat(R0, 1,length(tt));
    length(rr)
    length(R0_vector)
    length(tt)
    tt=tt(1:length(rr));
    length(tt)
    R0_vector=R0_vector(1:length(rr));
    length(R0_vector)

    T = table(transpose(rr), transpose(tt), transpose(R0_vector),'VariableNames',{'radius', 'time', 'R0'});
    writetable(T,sprintf('%s/radius_0.5_level_set_%s_epsilon_%s_alpha_%s.txt',indir, string(Nx),epsilon_name, alpha_name),'WriteMode','append')

    % M = [[rr; nan], [tt; nan], R0_vector]
end
