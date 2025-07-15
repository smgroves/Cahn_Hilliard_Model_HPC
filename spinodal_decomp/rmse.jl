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
# method = ARGS[4]
note= ARGS[5]
SLURM_ID=ARGS[6]

nx = GridSize
ny = nx

n_relax = 4
h = 1 / GridSize
m = 8
epsilon = m * h / (2 * sqrt(2) * atanh(0.9))

printphi = true
dt = 5.5e-6
max_it = 2000

if print
    dt_out = 10
else
    dt_out = 2000
end

id = @sprintf("SAV_MATLAB_%d_dt_%.2e_Nx_%d_%s_dt_out_%d%s", max_it, dt, GridSize, boundary, dt_out, note)
pathname1 = @sprintf("%s/out_MATLAB/%s_", outdir, id)

full_data1 = readdlm("$(pathname1)phi.csv", ',')
println(pathname1)

id = @sprintf("SAV_Python_%d_dt_%.2e_Nx_%d_%s_dtout_%d%s", max_it, dt, GridSize, boundary, dt_out, note)
pathname2 = @sprintf("%s/out_python/%s", outdir, id)
full_data2 = readdlm("$(pathname2)phi.csv", ',')
function parse_py_complex(x::String)
    x = strip(x, ['(', ')', ' '])         # Strip parens and whitespace
    if isempty(x)
        return NaN + 0im
    end
    x = replace(x, "j" => "im")
    try
        return Meta.parse(x) |> eval
    catch
        @warn "Could not parse complex string: $x"
        return NaN + 0im
    end
end

# Handle SubString and other types
parse_py_complex(x::SubString{String}) = parse_py_complex(String(x))
parse_py_complex(x::Float64) = ComplexF64(x)
parse_py_complex(x::Real) = ComplexF64(x)
parse_py_complex(x::Complex) = x
parse_py_complex(::Nothing) = NaN + 0im  # explicitly handle Nothing

full_data2 = real.(parse_py_complex.(full_data2))

println(pathname2)
#%%

#reshape full_data1 from [AxB,C] to [A,B,C]
timepoints = 201  # number of timepoints
tp1 = size(full_data1, 1) ÷ nx
tp2 = size(full_data2, 1) ÷ nx
if tp1 != tp2
    println("Timepoints do not match: $(tp1) vs $(tp2)")
end
full_data1 = reshape(full_data1, (nx, nx, tp1))
full_data1 = full_data1[:, :, 1:min(timepoints, tp2, tp1)]  # ensure we only take the first min(timepoints, tp1) timepoints
full_data2 = reshape(full_data2, (nx, nx, tp2))
full_data2 = full_data2[:, :, 1:min(timepoints, tp1, tp2)]  # ensure we only take the first min(timepoints, tp2) timepoints
# full_data2 = full_data2[:, :, 2:202]

# calculare L2 norm for each timepoint in phi
rmse = vec(sqrt.(mean((full_data1 - full_data2) .^ 2, dims=(1, 2))))
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

line = join(["$(title1)_$(title2)"; string.(rmse)], ",") * "\n"
open("$(outdir)/out_compare/rmse_SAV_MATLAB_Python.csv", "a") do io
    if !isfile("$(outdir)/out_python/rmse.csv") || filesize("$(outdir)/out_python/rmse.csv") == 0
        write(io, join(headers, ",") * "\n")
    end
    write(io, line)
end
