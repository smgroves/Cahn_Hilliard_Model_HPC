include("../CH_multigrid_solver_with_alpha_change_domain.jl")
using Dates
date_time = now()
nx = parse.(Int, ARGS[3])
ny = nx
tol = 1e-5
dt = 0.1 * (1 / 256)^2
epsilon = 0.0096

total_time = parse.(Float64, ARGS[4])
# total_time = 0.015
max_it = Int.(round(total_time / dt))
ns = 10
# CPC_width = 20
# experiments = .173 um = .108 of width
CPC_width = parse.(Float64, ARGS[1])
cohesin_width = parse.(Float64, ARGS[2])
println(CPC_width)
println(cohesin_width)
println("starting initialization")
phi = initialize_round_CPC(nx, nx, CPC_width=CPC_width, cohesin_width=cohesin_width)
# phi = initialize_geometric_CPC(nx, ny; CPC_width=20, cohesin_width=4)
outdir = "/project/g_bme-janeslab/sarah/julia_out/CPC_geometry/CPC_domain_0_2"
# outdir = "/scratch/xpz5km/Cahn_Hilliard_Model/CPC_geometry/CPC_domain_0_2"

println("starting ch solver")
println("cohesin=$(cohesin_width), CPC=$(CPC_width), epsilon=$(epsilon)")
phi = initialize_round_CPC_um(nx, nx, CPC_width=CPC_width, cohesin_width=cohesin_width, domain_width=6.4)
alpha = 0
time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_CPC_$(CPC_width)_cohesin_$(cohesin_width)_eps_$(epsilon)_alpha_$(alpha)_domain_0_2", check_dir=false, alpha=alpha)
open("/home/xpz5km/Cahn_Hilliard_Model/Job_specs.csv", "a", lock=false) do f
    writedlm(f, ["CPC_geometry_domain_0_2" "Julia" nx epsilon dt tol max_it 10000 time_passed], ",")
end



