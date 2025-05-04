include("../CH_multigrid_solver.jl")

nx = 1024
ny = 1024
tol = 1e-4
dt = 0.1 * (1 / 256)^2
m = 8
gam = m * (1 / 256) / (2 * sqrt(2) * atanh(0.9))
max_it_CH = 10000
# total_time = 0.03
# max_it = Int.(round(total_time / dt))
outdir = "/scratch/xpz5km/Cahn_Hilliard_Model/julia_out/CPC_geometry_1024/"
file = "$(outdir)/phi_1024_19660_0.0001__CPC_40_cohesin_16_eps_0.007504684956431058.txt"
num_lines = countlines(file)
println(num_lines)
#get initial from file
num_iterations_completed = Int.(round(num_lines/1024))
println(num_iterations_completed)
max_it = 19660 - num_iterations_completed
ns = 10
CPC_width = 40
# experiments = .173 um = .108 of width 
cohesin_width = 16
# phi = initialize_round_CPC(nx, nx, CPC_width=CPC_width, cohesin_width=cohesin_width)
# phi = initialize_geometric_CPC(nx, ny; CPC_width=20, cohesin_width=4)
phi = initialization_from_file("$(outdir)/tmp.txt", nx, ny, delim=' ') #from python code last_n_lines.py

time_passed = @elapsed main(phi, nx, tol, outdir, dt=dt, gam=gam, max_it=max_it, print_mass=false, print_e=false, print_r=true, overwrite=false, suffix="_CPC_$(CPC_width)_cohesin_$(cohesin_width)_eps_$(gam)", check_dir=false)
open("/home/xpz5km/Cahn_Hilliard_Model/Job_specs.csv", "a", lock=false) do f
    writedlm(f, ["CPC_geometry_partial" "Julia" nx gam dt tol max_it max_it_CH time_passed], ",")
end

