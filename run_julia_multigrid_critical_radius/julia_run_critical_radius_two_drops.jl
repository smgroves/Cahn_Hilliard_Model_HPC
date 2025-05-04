
using Dates
date_time = now()
include("../CH_multigrid_solver.jl")

function initialize_two_drops(nx, ny; h = 1/128, R0=0.1, gam=0.01)
    phi = zeros(Float64, nx, ny)
    x = h .* (0:nx-1)
    y = h .* (0:ny-1)
    xx, yy = meshgrid(x, y)
    y_center = 0.5
    x_center1 = 0.2
    x_center2 = .8

    R1 = @. sqrt((xx - x_center1)^2 + (yy - y_center)^2)
    R2 = @. sqrt((xx - x_center2)^2 + (yy - y_center)^2)

    phi = @.tanh((R0 .- R1) / (sqrt(2) * gam)) .+ @.tanh((R0 .- R2) / (sqrt(2) * gam)) .+ 1

    return phi
end

#ARGS: [R0, M, total_time]
nx = 512
tol = 1e-6
dt = 2.5e-5
R0 =parse.(Float64, ARGS[1])
M = parse.(Float64, ARGS[2]) # M = 8 is eps0 = 0.015
gam = M * (1 / 128) / (2 * sqrt(2) * atanh(0.9))
total_time = parse.(Int, ARGS[3])
max_it = Int.(total_time/dt)
max_it_CH = 10000
ns = 10
println(R0)
if R0 < 2.0
    println("starting initialization")
    phi = initialize_two_drops(nx, nx, h=1 / nx, R0 = R0, gam = gam)
    outdir="/project/g_bme-janeslab/sarah/julia_out/julia_out/critical_radius_two_drops"
    println("starting ch solver")

    time_passed = @elapsed main(phi, nx, tol, outdir, dt=dt, gam=gam, max_it=max_it, print_mass=false, print_e=false, overwrite=false, print_r=false, suffix="_R0_$(R0)_eps_$(round(gam, sigdigits = 5))_two_drops", check_dir=false)

else
    println("R0 was too big. Stopping run.")
end
