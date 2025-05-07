%source: https://github.com/nsbalbi/Spinodal-Decomposition
% adapted to use same initial conditions as our method 
% and to output energy, mass, and residual

function [phi_t,mass_t,E_t] = spinodal_decomp(D,gamma,options)
    % SPINODAL_DECOMP Generates and records Cahn-Hilliard spinodal
    % decomposition models using Euler's method.
    % 
    %    SPINODAL_DECOMP() generates a Cahn-Hilliad spindoal decomposition
    %    model with sample coefficients. The result is saved as an AVI named
    %    'spindoal_decomposition'.
    %
    %    SPINODAL_DECOMP(D,gamma) generates a Cahn-Hilliad spindoal 
    %    decomposition model with diffusion coefficient D and square length of
    %    transitional regions gamma. The result is saved as an AVI named 
    %    'spindoal_decomposition'. See "Inputs" for details.
    %    
    %    SPINODAL_DECOMP(...,PARAM1,VAL1,PARAM2,VAL2,...) are additional
    %    name-value pairs that can be used to change default function
    %    parameters. See "Parameters" for details.
    %
    % Inputs
    %    D: double. Diffusion coefficient. Default is 10
    %    gamma: double. Sqaure length of transitional regions between domains.
    %       Default is 5
    % 
    % Parameters
    %    'dt': double. Time increment between iterations of Euler's method. 
    %       (Default = 0.005)
    %    'GridSize': int. Edge length for the square grid used for the
    %       model. Note generation time increases exponentially with grid size.
    %       (Default = 200)
    %    'NumIterations': int. Total number of model iterations. Each iteration
    %       represents a step of Euler's method. (Default = 10000)
    %    'FrameSpacing': int. Number of iterations between captured frames, 
    %       only relevant if capture mode is standard. (Default = 10)
    %    'CaptureMode': char. Method of video cature. Possible inputs below.
    %         'standard' - Constant num of iterations between frames. (Default)
    %      'incremental' - Num iterations between frames increases over time
    %    'ImgStyle': char. Method of frame generation. Possible inputs below.
    %         'binary' - Concentrations are binarized to major species. (Default)
    %           'true' - Concentrations are mapped to the colormap.
    %    'Colormap': char. Colormap used for visualization. Supports all 
    %       default MATLAB colormaps (Default = 'pink')
    %    'FileName': char. Name of video (Default = 'spinodal_decomposition')
    %       
    % Examples
    %    Example 1: Generate and record default model.
    %       spinodal_decomp();
    %    Example 2: Generate and record model with custom constants.
    %       spinodal_decomp(20,3);
    %    Example 3: Generate and record model with custom constants and capture
    %    mode.
    %       spinodal_decomp(10,10,'CaptureMode','incremental');
    %    Example 4: Generate and record model with custom constants and
    %    multiple custom parameters.
    %       spinodal_decomp(10,10,...
    %                      'CaptureMode','incremental',...
    %                      'Colormap','jet',...
    %                      'ImgStyle','true',...
    %                      'NumIterations',25000);   

    arguments
    D double = 10
    gamma double = 5
    options.dt double = 0.005
    options.GridSize double = 200
    options.NumIterations double = 10000
    options.FrameSpacing double = 10
    options.CaptureMode char = 'standard'
    options.ImgStyle char = 'true'
    options.Colormap char = 'pink'
    options.FileName char = 'spinodal_decomposition'
    options.InputMatrix char = 'input.csv'
    options.InputType char = 'phi'
    options.ConstantColorbar = true
    options.write_phi = true
    options.write_residual = true
    options.ns = 1
    options.frame_stretch = 1 %default no frame stretching
    options.frame_per_sec = false % if true, each frame is 1 sec long * frame_stretch. If false, each frame is 1 frame length (30fps is default for MP4)
    end

    % Random Starting Concentrations (-1 and 1 represent different species)
    % u = 2*randi(2,options.GridSize)-3;
    u = readmatrix(options.InputMatrix);
    % if options.InputType == 'phi'
    %     u = (u+1)./2;
    % end


    writer = VideoWriter(options.FileName,'MPEG-4');
    writer.FrameRate = 30; % 30 frames per second

    open(writer);
    % figure(1)
    fig = figure('visible', 'off');

    % colormap(options.Colormap)
    colormap(redblue(100));
    image(u,'CDataMapping','scaled'); colorbar; axis square;
    set(gca,'FontSize',16);title('t = 0'); xlabel('x'); ylabel('y');
    colormap(interp1(1:100:1100,redbluecmap,1:1001)); %Expand redbluecmap to 1000 elements

    frame = getframe(fig);

    writeVideo(writer,frame);

    % Variables for incremental capture
    count = 1;
    frameStep = 1;

    % note that h2 is 1 because it is defined as ((xright-xleft)/nx)*((yright-yleft)/ny) and the domain is the same as gridsize (i.e. mesh size is 1)
    h2 = 1;
    E_t = zeros(options.NumIterations+1,1);
    mass_t = zeros(options.NumIterations+1,1);
    nx = options.GridSize;
    ny = nx;
    phi_t = zeros(nx,ny,options.NumIterations+1); %Initialize outputs
    if options.write_phi
        writematrix(u,sprintf('%s_phi.csv', options.FileName)); %write IC to file
    end

    for i = 1:options.NumIterations
        if rem(i, 1000) == 0
            i
        end
        curr_t=(i)*options.dt;
        mass_t(i) = sum(sum(u))/(h2*nx*ny);
        E_t(i) = discrete_energy(u,h2,nx,ny,gamma);

        u = iterate(u,D,gamma,options.dt);
        if options.write_phi
            if rem(i, options.ns) == 0
                writematrix(u,sprintf('%s_phi.csv', options.FileName),'WriteMode','append');
            end 
        end
        % Incremental video mode
        if (strcmp(options.CaptureMode,'incremental'))
            if (count == frameStep)
                if (strcmp(options.ImgStyle,'true'))
                    image(u,'CDataMapping','scaled');colorbar; axis square;
                else
                    uMod = round((u+1)/2); % Binarizes image
                    image(uMod,'CDataMapping','scaled');colorbar; axis square;
                end
                set(gca,'FontSize',16);title(['t = ',num2str(curr_t)]); xlabel('x'); ylabel('y');
                frame = getframe(fig);
                writeVideo(writer,frame);
                count = 0;
                frameStep = frameStep + 1;
            end
            count = count+1;
        % Standard video mode
        else
            if (mod(i,options.FrameSpacing) == 0)
                phi_t(:,:,i) = u;
                imAlpha=ones(size(u));
                imAlpha(isnan(u))=0;
                if (strcmp(options.ImgStyle,'true'))
                    image(u,'CDataMapping','scaled','AlphaData',imAlpha);colorbar; axis square;
                else
                    uMod = round((u+1)/2);
                    image(uMod,'CDataMapping','scaled','AlphaData',imAlpha);colorbar; axis square;
                end
                if options.ConstantColorbar
                    clim([-1, 1]);
                end
                colormap(interp1(1:100:1100,redbluecmap,1:1001)); %Expand redbluecmap to 1000 elements

                set(gca,'FontSize',16);title(['t = ',num2str(curr_t)]); xlabel('x'); ylabel('y');
                set(gca,'color',0*[1 1 1]);

                frame = getframe(fig);
                % Write this frame multiple times to stretch its duration
                if options.frame_per_sec
                    for repeat = 1:(options.frame_stretch * writer.FrameRate)
                        writeVideo(writer, frame);
                    end
                else
                    writeVideo(writer, frame);
                end

            end
        end
    end

    phi_t(:,:,options.NumIterations+1) = u;
    % final_phi = u;
    mass_t(options.NumIterations+1) = sum(sum(u))/(h2*nx*ny);
    E_t(options.NumIterations+1) = discrete_energy(u,h2,nx,ny,gamma);

    close(writer);

    %Normalize energy #note that we don't want to normalize mass because initial mass might be 0. 
    E_t = E_t/E_t(1);
    fprintf('Done!\n');

    % Forward Euler Method iteration of model
    function uOut = iterate(u,D,gamma,dt) 
        % Calculates laplacian of concentration field
        uLaplace = laplacian(u);
        % Calculates chemical potentials
        uMu = u.^3 - u - gamma*uLaplace;
        % Laplacian of chemical potentials
        muLaplace = laplacian(uMu);
        % Cahn-Hilliard Equation
        duT = D*muLaplace;
        % Foreward Euler Method
        uOut = u + dt*duT;
        if options.write_residual
            res2 = calc_residual(uOut,u,uMu,dt,D);
            fid = fopen(sprintf('%s_residual.txt',options.FileName), 'a+');
            fprintf(fid, '%f \n', res2);
            fclose(fid);
        end
    end 

end

%confirmed that residual is 0 for each iteration because this is an explicit method
function res2 = calc_residual(uNew, uOld, uMu, dt, D)
    s = size(uNew);
    Nx = s(1);
    rr = uMu; % slightly different than multigrid residual
    sor = D*laplacian(rr);
    rr = sor - (uNew - uOld)/dt;
    x = sum(sum(rr.*rr));
    res2 = sqrt(x / (Nx * Nx));
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

function E = discrete_energy(phi,hxy,gridx,gridy,eps2)
    %Local function for calculating the discrete energy across the domain
    f = @(x) 0.25*(x.^2-1).^2; %Define double-well free energy
    a = hxy*sum(sum(f(phi))); %Calculate chemical free energy
    sum_i = 0; %Initialize interfacial free energy in x
    for i = 1:gridx-1
        for j = 1:gridy
            sum_i = sum_i + (phi(i+1,j)-phi(i,j))^2;
        end
    end
    sum_j = 0; %Initialize interfacial free energy in y
    for i = 1:gridx
        for j = 1:gridy-1
            sum_j = sum_j + (phi(i,j+1)-phi(i,j))^2;
        end
    end
    E = a + 0.5*eps2*(sum_i+sum_j);
    end