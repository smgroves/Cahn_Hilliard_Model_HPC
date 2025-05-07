% indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/CPC_geometry/CPC_alpha_0";
% outdir = sprintf("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/%s",name);

% dt = 2.5e-5;
% dtout = 10;
function [] = kymograph_central_droplets_domain(indir, outdir, name, dt, dtout, transposed, cutoff, domain_um)

    mkdir(outdir) 

    phi = readmatrix(sprintf('%s/%s.txt', indir, name),'FileType','text');
    % phi = readmatrix("/Users/smgroves/Documents/GitHub/jlCHSolver/output.txt");
    phidims = size(phi);
    phidims(3) = phidims(1)/phidims(2); %Determine number of frames captured
    phidims(1) = phidims(2); %Determine size of square grid
    Nx = phidims(1);
    phi = reshape(phi,phidims(1),phidims(3),phidims(2)); %Reshape multidimensional array

    fig = figure('visible','off');

    x = (round(Nx/2));
    if cutoff == 0
        if transposed
            % if transposed, then the order is already correct for [x vs time]
            h=heatmap(phi(:,:,round(Nx/2)),'GridVisible','off');
        else
            % if not transposed, the order needs to change from [x, time, y] to [time, y, x] and then transpose time and y
            phi = shiftdim(phi,1); %Shift dimensions to move frames to the third dimension
            h=heatmap(transpose(phi(:,:,round(Nx/2))),'GridVisible','off');
        end
    else
        if transposed
            % if transposed, then the order is already correct for [x vs time]
            h=heatmap(phi(1:cutoff,:,round(Nx/2)),'GridVisible','off');
        else
            % if not transposed, the order needs to change from [x, time, y] to [time, y, x] and then transpose time and y
            phi = shiftdim(phi,1); %Shift dimensions to move frames to the third dimension
            h=heatmap(transpose(phi(1:cutoff,:,round(Nx/2))),'GridVisible','off');
        end
    end
    clim([-1, 1]);
    colormap(redblue(100));
    % colormap(h,parula);
    if cutoff == 0
        XLabels = 1:phidims(3);
    else
        XLabels = 1:cutoff;
    end
    CustomXLabels = string((XLabels-1)*dt*dtout);
    CustomXLabels(mod(XLabels,5000) ~= 0) = " ";
    % h.XDisplayLabels = CustomXLabels;
    % h.XLabel = "Time";

    YLabels = 1:phidims(2);
    % CustomYLabels = string(YLabels);

    % in um
    % 
    CustomYLabels = string(domain_um*(YLabels-x)/Nx);
    CustomYLabels(mod(domain_um*(YLabels-x)/Nx,.8) ~= 0) = " ";

    h.YDisplayLabels = CustomYLabels;
    h.YLabel = "Y (um)";
    h.Title = name;
    % set(gca, 'FontName', 'Arial');
    % set(gca,'YTickLabel',[]);
    % set(gca,'XTickLabel',[]);
    print(gcf,sprintf('%s/kymograph_x_%d_redblue_um.png', outdir, x),'-dpng')
    % print(gcf,sprintf('%s/kymograph_x_%d_redblue_um.pdf', outdir, x),'-dpdf','-vector')
end

function c = redblue(m)
    %   Adam Auton, 9th October 2009
    
    if nargin < 1, m = size(get(gcf,'colormap'),1); end
    
    if (mod(m,2) == 0)
        % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
        m1 = m*0.5;
        r = (0:m1-1)'/max(m1-1,1);
        g = r;
        r = [r; ones(m1,1)];
        g = [g; flipud(g)];
        b = flipud(r);
    else
        % From [0 0 1] to [1 1 1] to [1 0 0];
        m1 = floor(m*0.5);
        r = (0:m1-1)'/max(m1,1);
        g = r;
        r = [r; ones(m1+1,1)];
        g = [g; 1; flipud(g)];
        b = flipud(r);
    end
    
    c = [r g b]; 
end    