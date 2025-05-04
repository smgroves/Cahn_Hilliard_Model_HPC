indir = "/project/g_bme-janeslab/sarah/alpha_nonzero/CPC_geometry/CPC_alpha_-0.5";

f1 = figure;
hold on;
for CPC = ["0.12" "0.173" "0.22" "0.25" "0.3"]
    for cohesin = ["0.1" "0.2" "0.3"]
        maxes = [];
        CPC
        cohesin

        name =sprintf("phi_256_2000_1.0e-5__CPC_%s_cohesin_%s_eps_0.04_alpha_-0.5", CPC, cohesin);
        phi = readmatrix(sprintf('%s/%s.txt', indir, name),'FileType','text');
        phidims = size(phi);
        phidims(3) = phidims(1)/phidims(2); %Determine number of frames captured
        phidims(1) = phidims(2); %Determine size of square grid
        Nx = phidims(1);
        phi = reshape(phi,phidims(1),phidims(3),phidims(2)); %Reshape multidimensional array
        phi = shiftdim(phi,2); %Shift dimensions to move frames to the third dimension

        for t = 1:phidims(3)
            maxes = [maxes max(max(phi(:,:,t)))] 
        end
        eq_max = maxes(end);

        plot(1:phidims(3), maxes, '-', 'DisplayName', sprintf('Final Eq Max %f, CPC=%s, cohesin=%s', round(eq_max,3), CPC, cohesin));

    end
end

title('Maximum Phi (droplet phase) for alpha = -0.5');
legend show;
grid on;
hold off;

set(gcf, 'PaperSize', [11, 20])
orient(gcf,'landscape')
print(gcf,sprintf('%s/max.pdf', indir),"-dpdf",'-fillpage')


f2 = figure;
hold on;
for CPC = ["0.12" "0.173" "0.22" "0.25" "0.3"]
    for cohesin = ["0.1" "0.2" "0.3"]
        mins = [];
        CPC
        cohesin

        name =sprintf("phi_256_2000_1.0e-5__CPC_%s_cohesin_%s_eps_0.04_alpha_-0.5", CPC, cohesin);
        phi = readmatrix(sprintf('%s/%s.txt', indir, name),'FileType','text');
        phidims = size(phi);
        phidims(3) = phidims(1)/phidims(2); %Determine number of frames captured
        phidims(1) = phidims(2); %Determine size of square grid
        Nx = phidims(1);
        phi = reshape(phi,phidims(1),phidims(3),phidims(2)); %Reshape multidimensional array
        phi = shiftdim(phi,2); %Shift dimensions to move frames to the third dimension

        for t = 1:phidims(3)
            mins = [mins min(min(phi(:,:,t)))]
        end
        eq_min = mins(end);

        plot(1:phidims(3), mins, '-', 'DisplayName', sprintf('Final Eq Min %f, CPC=%s, cohesin=%s', round(eq_min,3), CPC, cohesin));

    end
end

title('Minimum Phi (droplet phase) for alpha = -0.5');
legend show;
grid on;
hold off;

set(gcf, 'PaperSize', [11, 20])
orient(gcf,'landscape')
print(gcf,sprintf('%s/min.pdf', indir),"-dpdf",'-fillpage')