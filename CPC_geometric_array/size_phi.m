% indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_0";
% outdir = sprintf("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/%s",name);

% dt = 2.5e-5;
% dtout = 10;
function [size_phi] = kymograph_central_droplets_domain(indir, name)

    phi = readmatrix(sprintf('%s/%s.txt', indir, name),'FileType','text');
    % phi = readmatrix("/Users/smgroves/Documents/GitHub/jlCHSolver/output.txt");
    phidims = size(phi);
    phidims(3) = phidims(1)/phidims(2); %Determine number of frames captured
    phidims(1) = phidims(2); %Determine size of square grid
    Nx = phidims(1);
    phi = reshape(phi,phidims(1),phidims(3),phidims(2)); %Reshape multidimensional array
    size_phi = size(phi)
end


