include("../CH_multigrid_solver.jl")
using Dates
date_time = now()

tol = 1e-5
dt = 0.1 * (1 / 256)^2
epsilon = 0.0075

CPC_width = parse.(Float64, ARGS[1])
cohesin_width = parse.(Float64, ARGS[2])
nx = parse.(Int, ARGS[3])
ny = nx
total_time = parse.(Float64, ARGS[4])
phi1 = parse.(Float64, ARGS[5])
phi2 = parse.(Float64, ARGS[6]) 
name = ARGS[7]


max_it = Int.(round(total_time / dt))
ns = 10

println(CPC_width)
println(cohesin_width)
println("starting initialization")
outdir = "/project/g_bme-janeslab/SarahG/julia_out/CPC_geometry/CPC_domain_0_1_e_0.0075_monse_quant"

println("starting ch solver")
println("cohesin=$(cohesin_width), CPC=$(CPC_width), epsilon=$(epsilon)")
phi = initialize_round_CPC_um(nx, nx, CPC_width=CPC_width, cohesin_width=cohesin_width, domain_width=3.2, phi1 = phi1, phi2 = phi2)
time_passed = @elapsed main(phi, nx, tol, outdir, dt=dt, gam=epsilon, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_CPC_$(CPC_width)_cohesin_$(cohesin_width)_eps_$(epsilon)_domain_0_1_$(name)", check_dir=false)

open("/home/xpz5km/Cahn_Hilliard_Model/Job_specs.csv", "a", lock=false) do f
    writedlm(f, ["CPC_geometry_e_0.0075_monse_quant" "Julia" nx epsilon dt tol max_it 10000 time_passed], ",")
end



