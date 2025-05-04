# using Pkg
# Pkg.add("DataFrames")
# Pkg.add("BenchmarkTools")
# Pkg.add("StaticArrays")
# Pkg.add("ProfileView")
# Pkg.add("DelimitedFiles")
# Pkg.add("LinearAlgebra")
# Pkg.add("Printf")
# Pkg.add("Dates")
using Dates
date_time = now()
include("./CH_multigrid_solver_with_alpha_v2.jl")
#ARGS: [R0, M, total_time, alpha,nx]
nx = parse.(Int, ARGS[5])
tol = 1e-6
dt = 2.5e-5
R0 =parse.(Float64, ARGS[1])
M = parse.(Float64, ARGS[2]) # M = 8 is eps0
gam = M * (1 / nx) / (2 * sqrt(2) * atanh(0.9))
total_time = parse.(Int, ARGS[3])
alpha=parse.(Float64,ARGS[4])
max_it = Int.(total_time/dt)
max_it_CH = 10000
ns = 10
println(R0)
println("starting initialization")
phi = initialization(nx, nx, method="droplet", h=1 / 128, R0 = R0, gam = gam)
outdir = "/scratch/xpz5km/Cahn_Hilliard_Model/julia_out/critical_radius_alpha"
println("starting ch solver")

time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=gam, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_R0_$(R0)_eps_$(round(gam, sigdigits = 5))_alpha_$(alpha)", check_dir=false, alpha=alpha)
open("/home/xpz5km/Cahn_Hilliard_Model/Job_specs.csv", "a", lock=false) do f
    writedlm(f, [date_time, "critical_radius_R0_$(R0)_eps_$(round(gam, sigdigits = 5))_alpha_$(alpha)" "Julia" nx gam dt tol max_it max_it_CH time_passed], ",")
end
