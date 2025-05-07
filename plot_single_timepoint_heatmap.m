function [] = plot_single_timepoint_heatmap(timestep, indir, name, suffix)
    %if you want a faster video, set frame_rate to something like 4

    % close all hidden
    % warning off
    const_colorbar = true;
    phi = readmatrix(sprintf('%s/%s.txt', indir, name),'FileType','text');
    % phi = readmatrix("/Users/smgroves/Documents/GitHub/jlCHSolver/output.txt");
    phidims = size(phi)
    phidims(3) = phidims(1)/phidims(2); %Determine number of frames captured
    phidims(1) = phidims(2); %Determine size of square grid
    phi = reshape(phi,phidims(1),phidims(3),phidims(2)); %Reshape multidimensional array
    phi = shiftdim(phi,2); %Shift dimensions to move frames to the third dimension
    myfig = figure();
    % hold on
    fig = figure('visible','off');
    h=heatmap(phi(:,:,timestep), 'CellLabelColor','none', 'GridVisible','off'); 


    XLabels = 1:phidims(2);
    % Convert each number in the array into a string
    CustomXLabels = string(XLabels);
    % Replace all but the fifth elements by spaces
    CustomXLabels(mod(XLabels,10) ~= 0) = " ";
    % Set the 'XDisplayLabels' property of the heatmap 
    % object 'h' to the custom x-axis tick labels
    h.XDisplayLabels = CustomXLabels;
    h.YDisplayLabels = CustomXLabels;

    colormap(redblue(100));
    % Set the color axis limits for the current frame
    set(gca,'FontSize',16);title(['t = ',num2str(timestep)]); xlabel('x'); ylabel('y');
    saveas(gca, sprintf('%s/%s%s.pdf', indir, name, suffix))
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


%{
module load matlab
indir="/project/g_bme-janeslab/SarahG/julia_out/CPC_geometry/CPC_domain_0_2_e_0.0075_noisy_cohesin/sd_0.11/individual_seeds"
name="phi_512_19661_1.0e-5__CPC_0.22_cohesin_0.09_eps_0.0075_alpha_0_domain_0_2"
suffix=""
matlab -nodisplay -r "plot_single_timepoint_heatmap(1, '$indir','$name', '$suffix');quit;"
%}