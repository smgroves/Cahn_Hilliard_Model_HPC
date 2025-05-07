function [t_out, phi_t, delta_mass_t, E_t] = CahnHilliard_SAV(phi0, varargin)
% This function uses the scalar auxiliary variable method to solve the 
% Cahn-Hilliard equation for a specified number of time steps of size dt.
% 
% INPUTS
    % phi0 = Initial field of chemical states in the domain, created by ch_initialization.
%
%NAME-VALUE PAIRS
    % t_iter = Number of time steps simulated. Default = 1e3.
    % dt = Time step. Default = 2.5e-5 characteristic times.
    %dt_out = Spacing of time steps output to phi_t as a multidimensional
    %   array (if less than 1e9 elements) or printed file (if greater than
    %   1e9 elements). Default = 1, to output every time step. 
    % m = Number of mesh points over which the interface exists. Default = 4.
    % epsilon2 = Squared interfacial transition distance; if specified,
    %   m will be overwritten. Default = nan (do not overwrite m).
    % boundary = Boundary conditions for the simulation:
    %   'periodic' (default) - flux on one domain border equals negative flux on the opposite border.
    %   'neumann' - zero flux on the domain borders.
    % domain = Vector of rightmost and leftmost grid points in x and y.
    %   Format: [xright xleft yright yleft]. Default = [1 0 1 0].
    % printphi = Logical to print phi to a file. Default = true.
    % pathname = Name of the path to which phi is printed. Default = 'cd'.
%
%OUTPUT
    % t_out = Time corresponding to the dt time step outputs.
    % phi_t = Multidimensional array of phi over t_out.
    % delta_mass_t = Vector of mass change over t_out.
    % E_t = Vector of total energy over t_out.
    
%% Set option defaults and parse inputs

    % Set parameter defaults
        default_t_iter = 1e3;
        default_dt = 2.5e-5;
        default_dt_out = 1;
        default_m = 4;
        default_epsilon2 = nan;
        default_boundary = 'periodic';
        default_domain = [1 0 1 0];
        default_printphi = false;
        default_pathname = 'cd';
        default_C0 = 0; %1
        default_Beta = 0;
        default_gamma0 = 0; %4

        % eta = 0.95
        % xi_flag = 1
        % Mob = 1
        
        CahnHilliard_SAV_parser = inputParser;

    % Set general criteria for inputs and name-value pairs
        valid_matrix = @(x) ismatrix(x);
        valid_integer = @(x) x-floor(x) == 0;
        valid_pos_num = @(x) isnumeric(x) && (x > 0);
        valid_boundary_type = @(x) strcmpi(x,'periodic') || strcmpi(x,'neumann');
        valid_domain_vector = @(x) length(x) == 4;
        valid_logical = @(x) islogical(x) || x == 1 || x == 0;
        valid_string = @(x) ischar(x) || isstring(x);

    % Set parser options and valid input criteria
        addRequired(CahnHilliard_SAV_parser,'phi0',valid_matrix);
        
        addParameter(CahnHilliard_SAV_parser,'t_iter',default_t_iter,valid_integer);
        addParameter(CahnHilliard_SAV_parser,'dt',default_dt,valid_pos_num);
        addParameter(CahnHilliard_SAV_parser,'dt_out',default_dt_out,valid_integer);
        addParameter(CahnHilliard_SAV_parser,'m',default_m,valid_integer);
        addParameter(CahnHilliard_SAV_parser,'epsilon2',default_epsilon2,valid_pos_num);
        addParameter(CahnHilliard_SAV_parser,'domain',default_domain,valid_domain_vector);
        addParameter(CahnHilliard_SAV_parser,'boundary',default_boundary,valid_boundary_type);
        addParameter(CahnHilliard_SAV_parser,'printphi',default_printphi,valid_logical);
        addParameter(CahnHilliard_SAV_parser,'pathname',default_pathname,valid_string);
        addParameter(CahnHilliard_SAV_parser,'C0',default_C0,valid_integer);
        addParameter(CahnHilliard_SAV_parser,'Beta',default_Beta,valid_integer);
        addParameter(CahnHilliard_SAV_parser,'gamma0',default_gamma0,valid_integer);

        parse(CahnHilliard_SAV_parser, phi0, varargin{:});
    
    % Extract parsed inputs
        phi0 = CahnHilliard_SAV_parser.Results.phi0;
        t_iter = CahnHilliard_SAV_parser.Results.t_iter;
        dt = CahnHilliard_SAV_parser.Results.dt;
        dt_out = CahnHilliard_SAV_parser.Results.dt_out;
        if dt_out > t_iter, error('Error: dt_out should not be greater than t_iter.'); end
        m = CahnHilliard_SAV_parser.Results.m;
        epsilon2 = CahnHilliard_SAV_parser.Results.epsilon2;
        boundary = CahnHilliard_SAV_parser.Results.boundary;
        xright = CahnHilliard_SAV_parser.Results.domain(1);
        xleft = CahnHilliard_SAV_parser.Results.domain(2);
        yright = CahnHilliard_SAV_parser.Results.domain(3);
        yleft = CahnHilliard_SAV_parser.Results.domain(4);
        printphi = CahnHilliard_SAV_parser.Results.printphi;
        pathname = CahnHilliard_SAV_parser.Results.pathname;
        C0 = CahnHilliard_SAV_parser.Results.C0;
        Beta = CahnHilliard_SAV_parser.Results.Beta;
        gamma0 = CahnHilliard_SAV_parser.Results.gamma0;

%% Define and initialize key simulation parameters

    [nx,ny] = size(phi0); % Define number of grid points in x and y
    Lx = xright-xleft; Ly = yright-yleft;

    % Decide on the solver's mesh spacing for NEUMANN vs PERIODIC
    %  - For Neumann: we will mirror the domain, so pass 2*hx and 2*hy into sav_solver.
    %  - For Periodic: keep as-is.
        if strcmpi(boundary,'neumann')
            Lx = 2*Lx;
            Ly = 2*Ly;
            nx = 2*nx;
            ny = 2*ny;
        elseif strcmpi(boundary,'periodic')
            Lx = Lx;
            Ly = Ly;
            nx = nx;
            ny = ny;
        end

    hx = Lx/nx; hy = Ly/ny;
    h2 = hx*hy; % Define mesh size

    if isnan(epsilon2)
        epsilon2 = h2*m^2/(2*sqrt(2)*atanh(0.9))^2; % Define ϵ^2 if not prespecified
    else
        m = sqrt((epsilon2*(2*sqrt(2)*atanh(0.9))^2)/h2); % Else overwrite m
        display(m);
    end

    k_x = 1i*[0:nx/2 -nx/2+1:-1]*(2*pi/Lx); k_y = 1i*[0:ny/2 -ny/2+1:-1]*(2*pi/Ly);
    k_xx = k_x.^2; k_yy = k_y.^2;
    [kxx,kyy] = meshgrid(k_xx,k_yy);
    k2 = kxx + kyy;
    k4 = k2.^2;

%% Initialization

    % Initialize chemical state and SAV state
        if strcmpi(boundary,'neumann')
            phi_old = ext(phi0); % Initialize chemical state with mirror extension
        elseif strcmpi(boundary,'periodic')
            phi_old = phi0; % Initialize chemical state
        end
        r_old = r0_fun(phi_old,hx,hy,C0); % Initialize sav state

    % Initialize output variables according to the output specifications
        n_timesteps = floor(t_iter/dt_out);
        downsampled = nx*ny*n_timesteps > 1e9; % Logical index for the need to downsample
        if printphi
            mass_t = zeros(n_timesteps+1,1);
            E_t = zeros(n_timesteps+1,1);
            if pathname == "cd"
                pathname = pwd;
            end
            Filename = strcat(pathname, 'phi.csv');
            % Note that this will overwrite the file if it already exists
            if strcmpi(boundary,'neumann')
                phi_old_out = extback(phi_old);
            elseif strcmpi(boundary,'periodic')
                phi_old_out = phi_old;
            end
            writematrix(phi_old_out, Filename, 'WriteMode', 'overwrite'); 
            phi_t = phi_old_out; %if printing out, just save the initial phi as phi_t so you don't get an error of ouput argument not assigned

        else
            if downsampled
                new_dt_out = ceil(nx*ny*t_iter/1e9); %we need to round up to ensure we have enough space
                fprintf("Variable phi_t is too large with dt_out = %4.0f. Downsampling to every %4.0f time steps\n", dt_out, new_dt_out)
                dt_out = new_dt_out;
                n_timesteps = floor(t_iter/dt_out);
            end
            if strcmpi(boundary,'neumann')
                phi_t = zeros(nx/2,ny/2,n_timesteps+1); 
                phi_old_out = extback(phi_old);
            elseif strcmpi(boundary,'periodic')
                phi_t = zeros(nx,ny,n_timesteps+1);
                phi_old_out = phi_old;
            end
            mass_t = zeros(n_timesteps+1,1);
            E_t = zeros(n_timesteps+1,1);
            phi_t(:,:,1) = phi_old_out;
        end
        mass_t(1) = ch_mass(phi_old_out,h2);
        E_t(1) = ch_discrete_energy(phi_old_out,h2,epsilon2);

%% Run SAV solver

    for i = 1:t_iter

        % Calculate current phi, r, mass and E
            [phi_new, r_new] = sav_solver(phi_old, r_old, ...
                hx, hy, k2, k4, dt, epsilon2, boundary, C0, Beta, gamma0);

        % Shrink the result back to the original domain size in phi_new_out for output
            if strcmpi(boundary,'neumann')
                phi_new_out = extback(phi_new);
            elseif strcmpi(boundary,'periodic')
                phi_new_out = phi_new;
            end

        % Calculate mass and energy according to the phi_new_out
            mass = ch_mass(phi_new_out,h2);
            E = ch_discrete_energy(phi_new_out,h2,epsilon2);

        % Store data according to the output specifications
            if mod(i,dt_out) == 0
                t_index = floor(i/dt_out)+1;
                if printphi
                % Write phi_new_out to file
                    writematrix(phi_new_out, Filename, 'WriteMode', 'append'); 
                else
                % Store as phi_t
                    phi_t(:,:,t_index) = phi_new_out;
                end
                % Store mass and energy
                mass_t(t_index) = mass;
                E_t(t_index) = E;
            end
            
        % Update iteration variables
            phi_old = phi_new;
            r_old = r_new;

        % Print percentage of completion
            if mod(i/t_iter*100,5) == 0
                fprintf('%3.0f percent complete\n',i/t_iter*100)
            end
    end

%% For post-processing

% Center mass and normalize energy to t == 0
    delta_mass_t = mass_t - mass_t(1);
    E_t = E_t/E_t(1);

% Output t_out vector for post-processing
    t_out = (0:1:n_timesteps)*dt*dt_out;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local functions

% Local function for calculating mass across the domain
    function mass = ch_mass(phi,h2)
        [nx,ny] = size(phi);
        mass = sum(sum(phi))/(h2*nx*ny);
    end

% Local function for "flip" extension
    function x_ext = ext(x)
        % Mirroring extension: takes x(nx, ny) -> x_ext(2*nx, 2*ny)
            [nx, ny] = size(x);
            x_ext = zeros(2*nx, 2*ny);
        
            % Original block
            x_ext(1:nx, 1:ny) = x;
        
            % Flip horizontally
            x_ext(1:nx, ny+1:2*ny) = x(:, end:-1:1);
        
            % Flip vertically
            x_ext(nx+1:2*nx, 1:ny) = x(end:-1:1, :);
        
            % Flip both
            x_ext(nx+1:2*nx, ny+1:2*ny) = x(end:-1:1, end:-1:1);
        end

% Local function for "flip" extension back
    function x_back = extback(x_ext)
        % Shrinks from 2*nx x 2*ny back to nx x ny (upper-left block)
            [nx_ext, ny_ext] = size(x_ext);
            nx = nx_ext/2;
            ny = ny_ext/2;
            x_back = x_ext(1:nx, 1:ny);
    end