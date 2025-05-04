
indir = "/project/g_bme-janeslab/SarahG/julia_out/critical_radius_updated_IC_256/";
everyR = 10000
f1 = figure;
hold on;
for R0 = ["0.117" "0.118"]
%for R0 = ["0.022","0.021","0.024","0.0215","0.0235","0.0205","0.02","0.01","0.025"]

    R0
    maxes = [];
    name =sprintf("phi_256_400000_1.0e-6__256_R0_%s_eps_0.015009", R0);
    phi = readmatrix(sprintf('%s/%s.txt', indir, name),'FileType','text');
    phidims = size(phi);
    phidims(3) = phidims(1)/phidims(2); %Determine number of frames captured
    phidims(1) = phidims(2); %Determine size of square grid
    Nx = phidims(1);
    phi = reshape(phi,phidims(1),phidims(3),phidims(2)); %Reshape multidimensional array
    phi = shiftdim(phi,2); %Shift dimensions to move frames to the third dimension

    for t = 1:everyR:phidims(3)
        maxes = [maxes max(max(phi(:,:,t)))]; 
    end
    eq_max = maxes(end);

    plot(1:everyR:phidims(3), maxes, '-', 'DisplayName', sprintf('Final Eq Max %f, R0 = %s', round(eq_max,3), R0));

end

title('Maximum Phi (droplet phase) for alpha = -0.5');
legend("Location", "southeast");
grid on;
hold off;

set(gcf, 'PaperSize', [8.5, 11])
orient(gcf,'landscape')
print(gcf,sprintf('%s/max_eps_0.015009_alpha_0.pdf', indir),"-dpdf",'-fillpage')


f2 = figure;
hold on;
for R0 = ["0.117" "0.118"]
%for R0 = ["0.07" "0.08" "0.09" "0.1" "0.11" "0.12" "0.13" "0.14"] %for e = 0.045
%for R0 = ["0.022","0.021","0.024","0.0215","0.0235","0.0205","0.02","0.01","0.025"]
    R0
    mins = [];

    name =sprintf("phi_256_400000_1.0e-6__256_R0_%s_eps_0.015009", R0);
    phi = readmatrix(sprintf('%s/%s.txt', indir, name),'FileType','text');
    phidims = size(phi);
    phidims(3) = phidims(1)/phidims(2); %Determine number of frames captured
    phidims(1) = phidims(2); %Determine size of square grid
    Nx = phidims(1);
    phi = reshape(phi,phidims(1),phidims(3),phidims(2)); %Reshape multidimensional array
    phi = shiftdim(phi,2); %Shift dimensions to move frames to the third dimension

    for t = 1:everyR:phidims(3)
        mins = [mins min(min(phi(:,:,t)))];
    end
    eq_min = mins(end);

    plot(1:everyR:phidims(3), mins, '-', 'DisplayName', sprintf('Final Eq Min %f, R0 = %s', round(eq_min,3), R0));

end

title('Minimum Phi (droplet phase) for alpha = -0.5');
%legend show;
legend("Location", "northeast");

grid on;
hold off;

set(gcf, 'PaperSize', [8.5, 11])
orient(gcf,'landscape')
print(gcf,sprintf('%s/min_eps_0.015009_alpha_0.pdf', indir),"-dpdf",'-fillpage')