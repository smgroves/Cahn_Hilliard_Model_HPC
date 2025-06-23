% indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/outputs/julia/figure_2_spinodal_decomp";
% name = "phi_32_10000_1.0e-6_";
% function for CPC sim videos with mins at the title
function [] = CHplotting_function_minutes(indir, name, dt, dtout, suffix, frame_rate)
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
            curr_t=(i-1)*dtout*dt*2.84*60;
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
            axis([0 size(phi,2) 0 size(phi,1)]);
            title(strcat('t = ',num2str(curr_t), ' min'));
            colormap(interp1(1:100:1100,redbluecmap,1:1001)); %Expand redbluecmap to 1000 elements
            colorbar;
            print(h, sprintf('%s/%s%s_tmp.png', indir, name, suffix), '-dpng', '-r300');
            % exportgraphics(h, 'frame.png', 'Resolution', 300);  % Save with high resolution
            img = imread(sprintf('%s/%s%s_tmp.png', indir, name, suffix));
            writeVideo(phi_movie, im2frame(img));
            % writeVideo(phi_movie,getframe(h));
            if mod(i,20)==0
                fprintf("%f%% \n",round(100*i/size(phi,3),2));
            end

        end
    end

    close(phi_movie);
end
