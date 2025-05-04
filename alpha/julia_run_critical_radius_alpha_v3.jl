using Dates
date_time = now()
include("./CH_multigrid_solver_with_alpha_v2.jl")
#ARGS: [R0, M, total_time, alpha,nx]
tol = 1e-6
dt = 2.5e-5
R0 =parse.(Float64, ARGS[1])
gam = parse.(Float64, ARGS[2])  # gam = M * (1 / nx) / (2 * sqrt(2) * atanh(0.9))
total_time = parse.(Float64, ARGS[3])
alpha=parse.(Float64,ARGS[4])
nx = parse.(Int, ARGS[5])

max_it = Int.(total_time/dt)
max_it_CH = 10000
ns = 10
println(R0)
println("starting initialization")
if alpha == -0.5
    beta = 0.679355124356581
    phi_1 = -1.2074996663997333
    phi_2 = 0.8741663330663997
elseif alpha == -0.1
    beta = 0.7024393586862704
    phi_1 = -1.0349986134211147
    phi_2 = 0.968331946754448
elseif alpha == 0.1
    beta = 0.7059312080025176
    phi_1 = -0.9683319467544478
    phi_2 =1.0349986134211147
elseif alpha == 0.5
    beta = 0.679366220486757
    phi_1 = -0.8741663330663997
    phi_2 = 1.2074996663997333
end

phi = initialization_from_function_v2(nx, nx, beta, phi_1, phi_2, R0=R0, gam = gam)
outdir = "/project/g_bme-janeslab/SarahG/julia_out/critical_radius_alpha_updated"
println("starting ch solver")

time_passed = @elapsed main_w_alpha(phi, nx, tol, outdir, dt=dt, gam=gam, max_it=max_it, print_mass=false, print_e=false, overwrite=false, suffix="_R0_$(R0)_eps_$(round(gam, sigdigits = 5))_alpha_$(alpha)", check_dir=false, alpha=alpha)
open("/home/xpz5km/Cahn_Hilliard_Model/Job_specs.csv", "a", lock=false) do f
    writedlm(f, [date_time, "critical_radius_R0_$(R0)_eps_$(round(gam, sigdigits = 5))_alpha_$(alpha)" "Julia" nx gam dt tol max_it max_it_CH time_passed], ",")
end
