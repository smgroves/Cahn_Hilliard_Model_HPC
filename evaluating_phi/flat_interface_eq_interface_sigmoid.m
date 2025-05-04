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
        writematrix(phi(:,:,phidims(3)),sprintf("%s/nx_%s_alpha_%s_eps_0.015009_equilibrium.csv", indir, nx, alpha))

    end
end
