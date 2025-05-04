
indir = "/project/g_bme-janeslab/SarahG/julia_out/critical_radius_alpha/flat_interface/";
everyR = 1000
f1 = figure;
hold on;
for nx = ["128", "256"]
    nx
    for alpha = ["-0.5", "-0.1", "0.1", "0.5"]
        alpha
        maxes = [];
        name =sprintf("phi_%s_200000_1.0e-6__nx_%s_eps_0.015009_alpha_%s", nx, nx, alpha);
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

        plot(1:everyR:phidims(3), maxes, '-', 'DisplayName', sprintf('Final Eq Max %f, Alpha = %s, nx = %s', round(eq_max,3), alpha, nx));
    end
end

title('Maximum Phi (droplet phase) for alpha = -0.5');
legend("Location", "southeast");
grid on;
hold off;

set(gcf, 'PaperSize', [8.5, 11])
orient(gcf,'landscape')
print(gcf,sprintf('%s/max_eps_0.015009_nonzero_alpha.pdf', indir),"-dpdf",'-fillpage')


f2 = figure;
hold on;
for nx = ["128", "256"]
    nx
    for alpha = ["-0.5", "-0.1", "0.1", "0.5"]
        alpha
        mins = [];

        name =sprintf("phi_%s_200000_1.0e-6__nx_%s_eps_0.015009_alpha_%s", nx, nx, alpha);
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

        plot(1:everyR:phidims(3), mins, '-', 'DisplayName', sprintf('Final Eq Min %f, Alpha = %s, nx = %s', round(eq_min,3), alpha, nx));
    end
end

title('Minimum Phi (droplet phase) for alpha = -0.5');
%legend show;
legend("Location", "northeast");

grid on;
hold off;

set(gcf, 'PaperSize', [8.5, 11])
orient(gcf,'landscape')
print(gcf,sprintf('%s/min_eps_0.015009_nonzero_alpha.pdf', indir),"-dpdf",'-fillpage')