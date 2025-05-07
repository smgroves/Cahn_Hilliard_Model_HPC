function [] = plot_IC_heatmap(indir, name, suffix, Nx, label_spacing)
    %if you want a faster video, set frame_rate to something like 4

    % close all hidden
    % warning off
    const_colorbar = true;
    % Open the file
    fid = fopen(sprintf('%s/%s.txt', indir, name), 'r');
    % Define the format specifier for 256 columns of floats
    formatSpec = repmat('%f', 1, Nx);  % Creates a format specifier with 256 '%f' for each column

    % Read the first n rows
    data = textscan(fid, formatSpec, Nx, 'Delimiter', ',', 'HeaderLines', 0);

    % Close the file
    fclose(fid);

    % Convert the cell array to a matrix for easier processing
    phi = cell2mat(data);
    myfig = figure();
    % hold on
    fig = figure('visible','off');
    h=heatmap(phi, 'CellLabelColor','none', 'GridVisible','off'); 

    phidims = size(phi);
    XLabels = 1:phidims(2);
    % Convert each number in the array into a string
    CustomXLabels = string(XLabels);
    % Replace all but the fifth elements by spaces
    CustomXLabels(mod(XLabels,label_spacing) ~= 0) = " ";
    % Set the 'XDisplayLabels' property of the heatmap 
    % object 'h' to the custom x-axis tick labels
    h.XDisplayLabels = CustomXLabels;
    h.YDisplayLabels = CustomXLabels;

    colormap(redblue(100));
    % Set the color axis limits for the current frame
    set(gca,'FontSize',16);title('Initial conditions'); xlabel('x'); ylabel('y');
    saveas(gca, sprintf('%s/%s%s_IC.pdf', indir, name, suffix))
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

