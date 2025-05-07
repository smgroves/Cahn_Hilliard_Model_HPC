include("../CH_multigrid_solver.jl")
using Dates
date_time = now()
nx = parse.(Int, ARGS[3])
ny = nx
tol = 1e-5
dt = 0.1 * (1 / 256)^2
# M = 8
# gam = M * (1 / 256) / (2 * sqrt(2) * atanh(0.9))
gam = 0.039117
max_it_CH = 10000
# total_time = 0.25
# max_it = Int.(round(total_time / dt))
max_it = 20000
ns = 10
# CPC_width = 20
# experiments = .173 um = .108 of width
CPC_width = parse.(Int, ARGS[1])
cohesin_width = parse.(Int, ARGS[2])
println(CPC_width)
println(cohesin_width)
println("starting initialization")
phi = initialize_round_CPC(nx, nx, CPC_width=CPC_width, cohesin_width=cohesin_width)
# phi = initialize_geometric_CPC(nx, ny; CPC_width=20, cohesin_width=4)
outdir = "/scratch/xpz5km/Cahn_Hilliard_Model/julia_out/critical_radius/CPC_geometry_eps_$(gam)"
println("starting ch solver")

time_passed = @elapsed main(phi, nx, tol, outdir, dt=dt, gam=gam, max_it=max_it, print_mass=false, print_e=false, print_r=false, overwrite=false, suffix="_CPC_$(CPC_width)_cohesin_$(cohesin_width)_eps_$(gam)", check_dir=false)
open("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/Job_specs.csv", "a", lock=false) do f
    writedlm(f, ["CPC_geometry_eps_$(gam)" "Julia" nx gam dt tol max_it max_it_CH time_passed], ",")
end

