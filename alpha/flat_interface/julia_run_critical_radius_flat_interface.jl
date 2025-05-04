using Dates
date_time = now()
include("./CH_multigrid_solver_with_alpha_v2.jl")
tol = 1e-6
dt = 2.5e-5

nx = parse.(Int, ARGS[1])
M = parse.(Float64, ARGS[2]) # M = 8 is eps0
gam = M * (1 / nx) / (2 * sqrt(2) * atanh(0.9))
total_time = parse.(Int, ARGS[3])
alpha=parse.(Float64,ARGS[4])

max_it = Int.(total_time/dt)
max_it_CH = 10000
ns = 10
println("starting initialization")
ny = nx

function initialization_tanh_transition(nx, ny, h; gam=0.1, interface_location=nx/2)
    phi = zeros(Float64, nx, ny)            # Initialize the phi array to zeros
    x = h .* (0:nx-1)                       # x-coordinates: scaled by step size h
    y = h .* (0:ny-1)                       # y-coordinates: scaled by step size h
    xx, yy = meshgrid(x, y)                 # Create a meshgrid of x and y
    
    # Create a transition along the x-axis
    phi = @.tanh((xx .- interface_location * h) / (sqrt(2) * gam))
    
    return phi
end

phi = initialization_tanh_transition(nx, nx, 1/256, gam = gam)

outdir = "/project/g_bme-janeslab/SarahG/julia_out/critical_radius_alpha/flat_interface"
println("starting ch solver")

time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=gam, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_nx_$(nx)_eps_$(round(gam, sigdigits = 5))_alpha_$(alpha)", check_dir=false, alpha=alpha)
open("/home/xpz5km/Cahn_Hilliard_Model/Job_specs.csv", "a", lock=false) do f
    writedlm(f, [date_time, "critical_radius_flat_interface_eps_$(round(gam, sigdigits = 5))_alpha_$(alpha)" "Julia" nx gam dt tol max_it max_it_CH time_passed], ",")
end

