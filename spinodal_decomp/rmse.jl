#%%
println("Starting...")
using CSV, DelimitedFiles
using Plots
using Printf
using Statistics
# import Pkg;
# Pkg.add("PyPlot");
# using PyPlot
# or use Makie for interactivity
# pyplot()
# --- Step 1: Read in the matrices ---
outdir="/project/g_bme-janeslab/SarahG/spinodal_decomp_06_2025"


GridSize = parse.(Int, ARGS[1])
boundary = ARGS[2]
print = lowercase(ARGS[3]) == "true"
method = ARGS[4]
note= ARGS[5]
SLURM_ID=ARGS[6]

nx = GridSize
ny = nx

n_relax = 4
h = 1 / GridSize
m = 8
epsilon = m * h / (2 * sqrt(2) * atanh(0.9))

printphi = true
dt = 5.5e-7
max_it = 20000
timepoints = 2001  # number of timepoints

if print
    dt_out = 10
else
    dt_out = 2000
end

#out_Julia/NMG_Julia_2000_dt_5.50e-06_Nx_512_neumann_dtout_1025p_phi
#out_MATLAB/NMG_MATLAB_2000_dt_5.50e-06_Nx_512_periodic_dt_out_1025p_phi
#out_Python/NMG_Python_20_dt_5.50e-06_Nx_512_neumann_dtout_150pphi
formatted_dt = @sprintf("%.2e", dt)

language = "Julia"

pathname1 = "$(outdir)/out_$(language)/NMG_$(language)_$(max_it)_dt_$(formatted_dt)_Nx_$(GridSize)_$(boundary)_dtout_$(dt_out)$(note)"
full_data1 = readdlm("$(pathname1)_phi.csv", ',')
println(pathname1)

# language = "Python"
pathname2 = "$(outdir)/out_$(language)/SAV_$(language)_$(max_it)_dt_$(formatted_dt)_Nx_$(GridSize)_$(boundary)_dtout_$(dt_out)$(note)"
full_data2 = readdlm("$(pathname2)_phi.csv", ',')
#for python code: there is an extra comma at the end that messes up the array size
# if all(x -> x == "" || x == missing, full_data2[:, end])
#     full_data2 = full_data2[:, 1:end-1]
# end

println(pathname2)
#%%

#reshape full_data1 from [AxB,C] to [A,B,C]
tp1 = size(full_data1, 1) ÷ nx
tp2 = size(full_data2, 1) ÷ nx
if tp1 != tp2
    println("Timepoints do not match: $(tp1) vs $(tp2)")
end


full_data1_reshaped = zeros(nx, nx, timepoints)
println(tp1)
println(tp2)
for i in 1:timepoints
    full_data1_reshaped[:, :, i] = full_data1[(nx*(i-1)+1):(nx*i), :]
end
full_data1_reshaped = full_data1_reshaped[:, :, 1:min(timepoints, tp2, tp1)]  # ensure we only take the first min(timepoints, tp1) timepoints

full_data2_reshaped = zeros(nx, nx, timepoints)
for i in 1:timepoints
    full_data2_reshaped[:, :, i] = full_data2[(nx*(i-1)+1):(nx*i), :]
end
full_data2_reshaped = full_data2_reshaped[:, :, 1:min(timepoints, tp1, tp2)]  # ensure we only take the first min(timepoints, tp2) timepoints

# calculare L2 norm for each timepoint in phi
rmse = vec(sqrt.(mean((full_data1_reshaped - full_data2_reshaped) .^ 2, dims=(1, 2))))
println(rmse)
ave_err = mean(rmse)

title1 = (basename(pathname1))
title2 = (basename(pathname2))
Plots.plot(rmse, label="RMSE", xlabel="Time Step", title="RMSE \n$(title1)\n$(title2)", titlefont=font(10), legend=:topleft, xlims=(1, min(timepoints, tp1, tp2)))
hline!([ave_err], linestyle=:dot, color=:black, label="Average RMSE = $(round(ave_err, digits = 4))")  # add horizontal dotted line at y = 0.0

#save FIGURE
println("saving figure")
savefig("$(outdir)/out_compare/plots/$(title1)_$(title2).pdf")

headers = ["compared_files";  string.(1:length(rmse))]
rmse_file = "$(outdir)/out_compare/rmse_NMG_SAV_MATLAB_dt_5.5e-7_v3.csv"
line = join(["$(title1)_$(title2)"; string.(rmse)], ",") * "\n"
open(rmse_file, "a") do io
    if !isfile(rmse_file) || filesize(rmse_file) == 0
        write(io, join(headers, ",") * "\n")
    end
    write(io, line)
end
