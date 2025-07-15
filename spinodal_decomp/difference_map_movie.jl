
using CSV, DelimitedFiles
using Plots
using Printf
using Statistics

default(colorbar=false)  # Turn off colorbars globally for cleaner animation

function difference_movie(difference, t_out, ny; dtframes=1, filename="difference_map", filetype="mp4", colorbar_type="default")
    """
    This function creates a red-white-blue video trajectory of chemical states using @animate.

    # Arguments
    - phi_file::String: File name for multidimensional array of chemical states.
    - t_out::Vector{Float64}: Vector of time steps corresponding to the third dimension of phi_t.
    - ny::Int: Number of mesh points in the y direction.

    # Optional Keyword Arguments
    - dtframes::Int: Number of dt time steps per frame. Default is 1. Note that this is a multiplier on dtout (i.e. if dtout = 10 and dtframes = 2, then every 20th timepoint is recorded in the movie). 
    - filename::String: Name of the movie file to be saved. Default is 'ch_movie'.
    - filetype::String: Format of the movie file ('mp4', 'gif'). Default is 'mp4'.
    - colorbar_type::String: Type of colorbar to be used ('default', 'initial_range'). Default is 'default'.

    # Output
    Saves a red-white-blue video.
    """

    # Validate file type
    valid_filetypes = ["mp4", "gif"]
    if !(filetype in valid_filetypes)
        error("Invalid filetype. Supported types are: mp4, gif")
    end

    # Initial color range
    initial_range = nothing

    anim = @animate for i in 1:dtframes:length(t_out)
        println(i)
        phi_temp = difference[:,:,i]

        # Determine color limits
        clims = if colorbar_type == "default"
            (-1, 1)
        elseif colorbar_type == "initial_range"
            initial_range === nothing && (initial_range = extrema(phi_temp))
            initial_range
        elseif colorbar_type == "variable"
        elseif colorbar_type == "variable"
            # extrema(phi_temp)
            vmin, vmax = extrema(phi_temp)
            vmax_abs = maximum(abs.((vmin, vmax)))
            if vmax_abs == 0
                vmax_abs = 1  # Avoid division by zero
            end
            (-vmax_abs, vmax_abs)
        else
            error("Invalid colorbar_type: $colorbar_type")
        end

        # Create the heatmap
        # title = @sprintf("%.3e", t_out[i])

        heatmap(phi_temp, y_flip=true, color=reverse(cgrad(:RdBu)), cbar=true, clims=clims, xlabel="", ylabel="", title="timestep = $(t_out[i])", xlim=(0, ny), ylim=(0, ny), aspect_ratio=:equal, dpi=300)
        if (i - 1) / (length(t_out) - 1) * 100 % 5 == 0
            println(@sprintf("%3.0f percent complete", (i - 1) / (length(t_out) - 1) * 100))
        end
    end #every dtframes <- redundant with above

    # Save animation
    if filetype == "mp4"
        mp4(anim, "$filename.mp4")
    elseif filetype == "gif"
        gif(anim, "$filename.gif")
    end

    println("Animation saved as $filename.$filetype")
end


println("Starting...")

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
dt = 5.5e-7
max_it = 20000

if print
    dt_out = 10
else
    dt_out = 2000
end

language = "Julia"
#NMG_MATLAB_2000_dt_5.50e-07_Nx_512_neumann_dt_out_1050p_phi
#NMG_Julia_20000_dt_5.50e-07_Nx_256_periodic_dtout_1075p_mass
id = @sprintf("SAV_%s_%d_dt_%.2e_Nx_%d_%s_dtout_%d%s", language,max_it, dt, GridSize, boundary, dt_out, note)
pathname1 = @sprintf("%s/out_julia/%s", outdir, id)

full_data1 = readdlm("$(pathname1)_phi.csv", ',')
println(pathname1)

id = @sprintf("NMG_%s_%d_dt_%.2e_Nx_%d_%s_dtout_%d%s", language, max_it, dt, GridSize, boundary, dt_out, note)
pathname2 = @sprintf("%s/out_juia/%s", outdir, id)
full_data2 = readdlm("$(pathname2)_phi.csv", ',')


#reshape full_data1 from [AxB,C] to [A,B,C]
timepoints = 2001  # number of timepoints
tp1 = size(full_data1, 1) ÷ nx
tp2 = size(full_data2, 1) ÷ nx
if tp1 != tp2
    println("Timepoints do not match: $(tp1) vs $(tp2)")
end
full_data1_reshaped = zeros(nx, nx, timepoints)

# full_data1 = reshape(full_data1, (nx, nx, tp1))
println(tp1)
println(tp2)
for i in 1:timepoints
    full_data1_reshaped[:, :, i] = full_data1[(nx*(i-1)+1):(nx*i), :]
end

full_data1_reshaped = full_data1_reshaped[:, :, 1:min(timepoints, tp2, tp1)]  # ensure we only take the first min(timepoints, tp1) timepoints
# full_data2 = reshape(full_data2, (nx, nx, tp2))
full_data2_reshaped = zeros(nx, nx, timepoints)

for i in 1:timepoints
    full_data2_reshaped[:, :, i] = full_data2[(nx*(i-1)+1):(nx*i), :]
end
full_data2_reshaped = full_data2_reshaped[:, :, 1:min(timepoints, tp1, tp2)]  # ensure we only take the first min(timepoints, tp2) timepoints

total_timesteps = min(timepoints, tp1, tp2)

title1 = (basename(pathname1))
title2 = (basename(pathname2))
println(size((full_data1_reshaped - full_data2_reshaped)))
difference_movie((full_data1_reshaped - full_data2_reshaped), 1:total_timesteps, nx; dtframes=1, filename="$(outdir)/out_julia/plots/v2_$(title1)_$(title2)", filetype="mp4", colorbar_type="default")