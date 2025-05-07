% indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/outputs/julia/figure_2_spinodal_decomp";
% name = "phi_32_10000_1.0e-6_";

function [] = CHplotting_function(indir, name, dt, dtout, suffix, frame_rate)
    %if you want a faster video, set frame_rate to something like 4

    % close all hidden
    % warning off
    const_colorbar = true;
    plot_type = "heatmap";
    phi = readmatrix(sprintf('%s/%s.txt', indir, name),'FileType','text');
    % phi = readmatrix("/Users/smgroves/Documents/GitHub/jlCHSolver/output.txt");
    phidims = size(phi)
    phidims(3) = phidims(1)/phidims(2); %Determine number of frames captured
    phidims(1) = phidims(2); %Determine size of square grid
    phi = reshape(phi,phidims(1),phidims(3),phidims(2)); %Reshape multidimensional array
    phi = shiftdim(phi,2); %Shift dimensions to move frames to the third dimension
    myfig = figure();
    % hold on
    if const_colorbar == true
        %note that Rivanna doesn't have mp4 processor
        phi_movie = VideoWriter(sprintf('%s/%s%s.avi', indir, name, suffix),'Motion JPEG AVI');
    else
        phi_movie = VideoWriter(sprintf('%s/_%s%s.avi', indir, name, suffix),'Motion JPEG AVI');
    end
    phi_movie.Quality = 100; %Highest quality video
    open(phi_movie);
    % Set consistent color axis limits for the entire movie
    clim([-1, 1]);
    % for i = 1:10
    for i = 1:phidims(3)
        if mod(i-1,frame_rate) == 0
            curr_t=(i-1)*dtout*dt;
            % % phi_tmp = circshift(phi_tmp,int64(phidims(1)/2),1); %%%%%%%%ADDED TO CHECK OUTPUT AFTER CIRCSHIFT
            h = figure('visible','off');
            image(phi(:,:,i),'CDataMapping','scaled');
            colorbar; 
            axis square;
            % if colorbar_type == "default"
            clim([-1, 1]);
            % elseif colorbar_type == "initial_range"
            %     clim([min(phi(:,:,1),[],'all') max(phi(:,:,1),[],'all')]); %Set color axis to initial conditions
            % end
            g = gca;
            %Scale and center display roughly to the size of the mesh,
            %with extra space horizontally for the color bar
            % g.Position = [0.5-size(phi,2)/1000 0.5-size(phi,1)/1000 ...
            %     2.2*size(phi,2)/1000 2*size(phi,1)/1000];
            axis([0 size(phi,2) 0 size(phi,1)]);
            title(strcat('t = ',num2str(curr_t)));
            colormap(interp1(1:100:1100,redbluecmap,1:1001)); %Expand redbluecmap to 1000 elements
            colorbar;
            writeVideo(phi_movie,getframe(h));
            if mod(i,20)==0
                fprintf("%f%% \n",round(100*i/size(phi,3),2));
            end

            % fig = figure('visible','off');
            % axis square;
            % if plot_type == "surf"
            %     surf(phimp,'EdgeColor','none');colorbar;
            %     view(2);
        
            % elseif plot_type == "contourf"
            %     contourf(phi_tmp); colorbar; axis square;
            %     view(2);
        
            % elseif plot_type == "surf3d"
            %     surfc(phi_tmp); colorbar; axis square;
            %     ylim([0, 2])
            % elseif plot_type == "surfside"
            %     surfc(phi_tmp); colorbar; axis square;
            %     view(90,0)
            %     ylim([0,2])
            % elseif plot_type == "heatmap"
            %     h=heatmap(phi_tmp, 'CellLabelColor','none', 'GridVisible','off'); 
            % end

            % % contour(:,:,i); colorbar; axis square;
            % if plot_type == "heatmap"
            %     XLabels = 1:phidims(2);
            %     % Convert each number in the array into a string
            %     CustomXLabels = string(XLabels);
            %     % Replace all but the fifth elements by spaces
            %     CustomXLabels(mod(XLabels,10) ~= 0) = " ";
            %     % Set the 'XDisplayLabels' property of the heatmap 
            %     % object 'h' to the custom x-axis tick labels
            %     h.XDisplayLabels = CustomXLabels;
            %     h.YDisplayLabels = CustomXLabels;

            % else
            %     axis([1 size(phi,1) 1 size(phi,2)]);
            % end
            % % colormap(redblue(100));
            % colormap(interp1(1:100:1100,redbluecmap,1:1001)); %Expand redbluecmap to 1000 elements

            % % Set the color axis limits for the current frame
            % set(gca,'FontSize',16);title(['t = ',num2str(curr_t)]); xlabel('x'); ylabel('y');
            % if const_colorbar == true
            %     clim([-1, 1]);
            % end
            % frame = getframe(fig); % Capture the current figure window as a frame
            % writeVideo(phi_movie, frame); % Add the frame to the video
            % % if mod(i,20)==0
            % %     close all;
            % % end
            % if mod(i,20*frame_rate)==0
            %     fprintf("%f%% \n",round(100*i/phidims(3),2));
            % end
        end
    end

    close(phi_movie);
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