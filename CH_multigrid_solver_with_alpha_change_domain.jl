# Julia implementation of a Multigrid solver for the Cahn-Hilliard equation
# Author: Sarah Groves
# June 7, 2024

using DataFrames
using DelimitedFiles
using LinearAlgebra
using BenchmarkTools
using StaticArrays
using Printf
import Random
import Pkg; Pkg.add("Distributions")
using Random, Distributions

Random.seed!(1234) #note that when using a random seed, you must RESTART the REPL or run in a new instance for the runs to be the same. 

const c_relax::Int = 2  # number of SMOOTH relaxation operations defined as a global variable
const xleft = 0.0  # left x-coordinate defined as a global variable
const xright = 2.0  # right x-coordinate defined as a global variable
const yleft = 0.0  # left y-coordinate defined as a global variable
const yright = 2.0  # right y-coordinate defined as a global variable

# laplacian function: laplacian(m, nx, ny, h2)
function laplace(a, nxt, nyt)

    lap_a = zeros(Float64, nxt, nyt)
    ht2 = ((xright - xleft) / nxt)^2
    for i in 1:nxt
        for j in 1:nyt
            if i > 1
                dadx_L = (a[i, j] - a[i-1, j])
            else
                dadx_L = 0
            end
            if i < nxt
                dadx_R = (a[i+1, j] - a[i, j])
            else
                dadx_R = 0
            end
            if j > 1
                dady_B = (a[i, j] - a[i, j-1])
            else
                dady_B = 0
            end
            if j < nyt
                dady_T = (a[i, j+1] - a[i, j])
            else
                dady_T = 0
            end
            lap_a[i, j] = (dadx_R - dadx_L + dady_T - dady_B) / ht2
        end
    end
    return lap_a
end

function print_mat(file, matrix)
    open(file, "a", lock=false) do f
        for i = axes(matrix, 1)
            for j = axes(matrix, 2)
                Printf.@printf(f, "%16.15f ", matrix[i, j])
            end
            println(f)
        end
    end
end

function source(c_old, nx, ny, dt)
    src_mu = zeros(Float64, nx, ny)
    src_c = zeros(Float64, nx, ny)
    ct = laplace(c_old, nx, ny)
    for i in 1:nx
        for j in 1:ny
            src_c[i, j] = c_old[i, j] / dt - ct[i, j]
            src_mu[i, j] = 0
        end
    end
    return src_c, src_mu
end

function relax(c_new, mu_new, su, sw, nxt, nyt, c_relax, xright, xleft, dt, Cahn, alpha)
    ht2 = ((xright - xleft) / nxt)^2
    # a = MVector{4,Float64}(undef)
    # f = MVector{2,Float64}(undef)
    a = zeros(Float64, 4)
    f = zeros(Float64, 2)
    for iter in 1:c_relax
        for i in 1:nxt
            for j in 1:nyt
                if i > 1 && i < nxt
                    x_fac = 2.0
                else
                    x_fac = 1.0
                end
                if j > 1 && j < nyt
                    y_fac = 2.0
                else
                    y_fac = 1.0
                end
                a[1] = 1 / dt
                a[2] = (x_fac + y_fac) / ht2
                #a[2] = -(x_fac + y_fac) * Cahn / ht2 - d2f(c_new[i][j]);
                a[3] = -(x_fac + y_fac) * Cahn / ht2 - 3 * (c_new[i, j])^2
                cnew_val = (c_new[i, j])
                d2f = -3 * (c_new[i, j])^2

                a[4] = 1.0

                f[1] = su[i, j]
                f[2] = sw[i, j] - 2 * (c_new[i, j])^3 - alpha * (c_new[i, j])^2 + alpha ## UPDATED

                if i > 1
                    f[1] += mu_new[i-1, j] / ht2
                    f[2] -= Cahn * c_new[i-1, j] / ht2
                end
                if i < nxt
                    f[1] += mu_new[i+1, j] / ht2
                    f[2] -= Cahn * c_new[i+1, j] / ht2
                end
                if j > 1
                    f[1] += mu_new[i, j-1] / ht2
                    f[2] -= Cahn * c_new[i, j-1] / ht2
                end
                if j < nyt
                    f[1] += mu_new[i, j+1] / ht2
                    f[2] -= Cahn * c_new[i, j+1] / ht2
                end
                det = a[1] * a[4] - a[2] * a[3]
                c_new[i, j] = (a[4] * f[1] - a[2] * f[2]) / det
                mu_new[i, j] = (-a[3] * f[1] + a[1] * f[2]) / det

            end
        end
    end
    return c_new, mu_new
end

function restrict_ch(uf, vf, nxc, nyc)
    uc = zeros(Float64, round(Int64, nxc), round(Int64, nyc))
    vc = zeros(Float64, round(Int64, nxc), round(Int64, nyc))
    for i in 1:nxc
        for j in 1:nyc
            uc[i, j] = 0.25 * (uf[round(Int, 2 * i - 1), round(Int, 2 * j - 1)] + uf[round(Int, 2 * i - 1), round(Int, 2 * j)] + uf[round(Int, 2 * i), round(Int, 2 * j - 1)] + uf[round(Int, 2 * i), round(Int, 2 * j)])
            vc[i, j] = 0.25 * (vf[round(Int, 2 * i - 1), round(Int, 2 * j - 1)] + vf[round(Int, 2 * i - 1), round(Int, 2 * j)] + vf[round(Int, 2 * i), round(Int, 2 * j - 1)] + vf[round(Int, 2 * i), round(Int, 2 * j)])
        end
    end
    return uc, vc
end

function nonL(c_new, mu_new, nxt, nyt, dt, Cahn, alpha)
    ru = zeros(Float64, nxt, nyt)
    rw = zeros(Float64, nxt, nyt)
    lap_c = laplace(c_new, nxt, nyt)
    lap_mu = laplace(mu_new, nxt, nyt)
    for i in 1:nxt
        for j in 1:nyt
            ru[i, j] = c_new[i, j] / dt - lap_mu[i, j]
            rw[i, j] = mu_new[i, j] / dt - c_new[i, j]^3 + Cahn * lap_c[i, j] + alpha * (c_new[i, j])^2 - alpha
        end
    end
    return ru, rw
end

# df(c) function: c.^3
# d2f(c) function: 3*c.^2
function defect(uf_new, wf_new, suf, swf, nxf, nyf, uc_new, wc_new, nxc, nyc, dt, Cahn, alpha)
    ruc, rwc = nonL(uc_new, wc_new, nxc, nyc, dt, Cahn, alpha)
    ruf, rwf = nonL(uf_new, wf_new, nxf, nyf, dt, Cahn, alpha)
    ruf = suf - ruf
    rwf = swf - rwf
    rruf, rrwf = restrict_ch(ruf, rwf, nxc, nyc)
    duc = ruc + rruf
    dwc = rwc + rrwf
    return duc, dwc
end

function prolong_ch(uc, vc, nxc, nyc)
    uf = zeros(Float64, 2 * nxc, 2 * nyc)
    vf = zeros(Float64, 2 * nxc, 2 * nyc)
    for i in 1:nxc
        for j in 1:nyc
            uf[2*i-1, 2*j-1] = uc[i, j]
            uf[2*i-1, 2*j] = uc[i, j]
            uf[2*i, 2*j-1] = uc[i, j]
            uf[2*i, 2*j] = uc[i, j]
            vf[2*i-1, 2*j-1] = vc[i, j]
            vf[2*i-1, 2*j] = vc[i, j]
            vf[2*i, 2*j-1] = vc[i, j]
            vf[2*i, 2*j] = vc[i, j]
        end
    end
    return uf, vf
end

function vcycle(uf_new, wf_new, su, sw, nxf, nyf, ilevel, c_relax, xright, xleft, dt, Cahn, n_level, alpha)
    uf_new, wf_new = relax(uf_new, wf_new, su, sw, nxf, nyf, c_relax, xright,
        xleft, dt, Cahn, alpha)
    if ilevel < n_level
        nxc = trunc(Int64, nxf / 2)
        nyc = trunc(Int64, nyf / 2)
        uc_new, wc_new = restrict_ch(uf_new, wf_new, nxc, nyc)
        duc, dwc = defect(uf_new, wf_new, su, sw, nxf, nyf, uc_new, wc_new, nxc, nyc, dt, Cahn, alpha)

        uc_def = copy(uc_new)
        wc_def = copy(wc_new)

        uc_def, wc_def = vcycle(uc_def, wc_def, duc, dwc, nxc, nyc, ilevel + 1, c_relax, xright, xleft, dt, Cahn, n_level, alpha)


        uc_def = uc_def - uc_new
        wc_def = wc_def - wc_new

        uf_def, wf_def = prolong_ch(uc_def, wc_def, nxc, nyc)

        uf_new = uf_new + uf_def
        wf_new = wf_new + wf_def

        uf_new, wf_new = relax(uf_new, wf_new, su, sw, nxf, nyf, c_relax, xright,
            xleft, dt, Cahn, alpha)

    end
    return uf_new, wf_new
end

function error2(c_old, c_new, mu, nxt, nyt, dt)
    rr = zeros(Float64, nxt, nyt)
    x = 0.0
    for i in 1:nxt
        for j in 1:nyt
            rr[i, j] = mu[i, j] - c_old[i, j]
        end
    end

    sor = laplace(rr, nxt, nyt)
    for i in 1:nxt
        for j in 1:nyt
            rr[i, j] = sor[i, j] - (c_new[i, j] - c_old[i, j]) / dt
        end
    end

    for i in 1:nxt
        for j in 1:nyt
            x += rr[i, j]^2
        end
    end
    res2 = sqrt(x / (nxt * nyt))
    return res2
end

function initialize_geometric_CPC(nx, ny, CPC_width=20, cohesin_width=4)
    phi = zeros(Float64, nx, ny)
    # CPC_width = 20
    # cohesin_width = 4
    for i in 1:nx
        for j in 1:ny
            if i > round(nx / 2) - CPC_width && i < round(nx / 2) + CPC_width
                if j > round(ny / 2) - CPC_width && j < round(ny / 2) + CPC_width
                    phi[i, j] = 1
                elseif i > round(nx / 2) - cohesin_width && i < round(nx / 2) + cohesin_width
                    phi[i, j] = 1
                else
                    phi[i, j] = -1
                end
            else
                phi[i, j] = -1
            end
        end
    end
    return phi
end

function initialize_round_CPC_um(nx, ny; CPC_width=0.173, cohesin_width=0.2, domain_width=3.2, c1 = -1.0, c2 = 1.0)
    CPC_radius_grid_points = nx * CPC_width / domain_width
    cohesin_width_grid_points = nx * cohesin_width / domain_width
    cohesin_half_width = cohesin_width_grid_points / 2
    # Create an empty matrix filled with -1
    phi = fill(c1, nx, ny)

    # Define the center of the matrix
    center = (nx) / 2

    # Loop through each element of the matrix
    for i in 1:nx
        for j in 1:ny
            # Calculate the distance from the center
            distance = norm([i - center, j - center])

            # Check if the distance is less than or equal to CPC_width
            if distance <= CPC_radius_grid_points
                phi[i, j] = c2
            elseif i > round((nx) / 2) - cohesin_half_width && i < round((nx) / 2) + cohesin_half_width
                phi[i, j] = c2
            end
        end
    end
    return phi
end

function initialize_round_CPC_noisy_cohesin_um(nx, ny; CPC_width=0.173, cohesin_width=0.2, domain_width=3.2, sd = 0.5, seed = 1234)
    Random.seed!(seed) 
    CPC_radius_grid_points = nx * CPC_width / domain_width
    # Create an empty matrix filled with -1
    phi = fill(-1.0, nx, ny)

    # Define the center of the matrix
    center = (nx) / 2

    # Loop through each element of the matrix
    for j in 1:ny
        noisy_cohesin_width = (rand(Normal(cohesin_width, cohesin_width*sd),1)[1])
        if noisy_cohesin_width < 0
            noisy_cohesin_width = 0
        end
        cohesin_width_grid_points = nx * noisy_cohesin_width / domain_width
        cohesin_half_width = cohesin_width_grid_points / 2

        for i in 1:nx
            # Calculate the distance from the center
            distance = norm([i - center, j - center])

            # Check if the distance is less than or equal to CPC_width
            if distance <= CPC_radius_grid_points
                phi[i, j] = 1.0
            elseif i > round((nx) / 2) - cohesin_half_width && i < round((nx) / 2) + cohesin_half_width
                phi[i, j] = 1.0
            end
        end
    end
    return phi
end

function initialize_round_CPC(nx, ny; CPC_width=10, cohesin_width=4)
    # Create an empty matrix filled with -1
    phi = fill(-1.0, nx, ny)

    # Define the center of the matrix
    center = (nx) / 2

    # Loop through each element of the matrix
    for i in 1:nx
        for j in 1:ny
            # Calculate the distance from the center
            distance = norm([i - center, j - center])

            # Check if the distance is less than or equal to CPC_width
            if distance <= CPC_width
                phi[i, j] = 1.0
            elseif i > round((nx) / 2) - cohesin_width && i < round((nx) / 2) + cohesin_width
                phi[i, j] = 1.0
            end
        end
    end
    return phi
end

function meshgrid(x, y)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

function initialization_from_function(nx, ny, h; R0=0.1, gam=0.01)
    phi = zeros(Float64, nx, ny)
    x = h .* (0:nx-1)
    y = h .* (0:ny-1)
    xx, yy = meshgrid(x, y)
    R = @.sqrt((xx - 0.5)^2 + (yy - 0.5)^2)
    # eps_c = gam
    # delta = eps_c * sqrt(2)
    # psi0 = 0.5 * (1 .+ @.tanh((R0 .- R) / (2 * delta)))
    # phi = 2 .* psi0 .- 1    # psi0=(phi0+1)/2
    phi = @.tanh((R0 .- R) / (sqrt(2) * gam))

    return phi
end

function initialization_random(nx, ny)
    return (1 .- 2 .* rand(nx, ny))
end

function initialization_spinodal(nx, ny)
    return (rand([-1.0, 1.0], nx, ny))
end

function initialization_from_file(file, nx, ny; delim=',', transpose_matrix=false)
    phi = readdlm(file, delim, Float64)
    if size(phi) != (nx, ny)
        print("Warning: phi from file is wrong size: $(size(phi)) Expected: $(nx), $(ny)")
    end
    if transpose_matrix
        phi = transpose(phi)
    end
    return phi
end

function initialization(nx, ny; method="spinodal", initial_file="", delim=",", h=1 / 128, R0=0.1, gam=0.01)
    if method == "random"
        oc = initialization_random(nx, ny)
    elseif method == "droplet"
        oc = initialization_from_function(nx, ny, h, R0=R0, gam=gam)
    elseif method == "geometric"
        oc = initialize_geometric_CPC(nx, ny)
    elseif method == "file"
        oc = initialization_from_file(initial_file, nx, ny, delim=delim)
    elseif method == "spinodal"
        oc = initialization_spinodal(nx, ny)
    else
        println("Warning: initialize must be one of [random, droplet, geometric, file, spinodal].")
    end
    return oc
end

function cahn(c_old, c_new, mu, nx, ny, dt, max_it_CH, tol, c_relax, xright, xleft, Cahn, n_level, alpha; suffix="")
    it_mg2 = 0
    resid2 = 1
    sc, smu = source(c_old, nx, ny, dt)

    while resid2 > tol && it_mg2 < max_it_CH

        c_new, mu = vcycle(c_new, mu, sc, smu, nx, ny, 1, c_relax, xright, xleft, dt, Cahn, n_level, alpha)

        resid2 = error2(c_old, c_new, mu, nx, ny, dt)
        print_mat("$(outdir)/residual_$(nx)_$(tol)_$(suffix).txt", resid2)
        it_mg2 += 1
    end

    return c_new
end

#will need to ensure matrices are floats not ints
function convert_matrix(matrix)
    if eltype(matrix) == Float64
        return matrix  # Matrix is already of type Float64, no need to convert
    elseif eltype(matrix) == Int
        return convert(Matrix{Float64}, matrix)  # Convert to Matrix{Float64}
    else
        throw(ArgumentError("Unsupported matrix element type"))
    end
end

function calculate_mass(phi, h2, nx, ny)
    ave_mass = sum(phi) / (h2 * nx * ny)
    return ave_mass
end

function f(phi)
    fphi = (1 / 4) .* ((phi .^ 2) .- 1) .^ 2
    return fphi
end

function calculate_discrete_energy(phi, h2, nx, ny, Cahn)
    a = h2 * sum(f(phi))
    sum_i = 0.0
    for i in 1:(nx-1)
        for j in 1:ny
            sum_i += (phi[i+1, j] - phi[i, j])^2

        end
    end
    b = (Cahn / 2) * sum_i
    sum_j = 0.0
    for i in 1:nx
        for j in 1:(ny-1)
            sum_j += (phi[i, j+1] - phi[i, j])^2
        end
    end
    c = (Cahn / 2) * sum_j
    E = a + b + c
    return E
end

function calculate_discrete_norm_energy(phi, phi0, h2, nx, ny, Cahn)
    E0 = calculate_discrete_energy(phi0, h2, nx, ny, Cahn)
    E = calculate_discrete_energy(phi, h2, nx, ny, Cahn)
    return E / E0
end

function main_w_alpha(oc, nx, tol, outdir; max_it=1000, max_it_CH=10000, suffix="", overwrite=true, print_phi=true, print_mass=false, print_e=false, dt=2.5e-5, M=8, ns=10, gam=0.0, alpha=0.0, check_dir=true)
    while true
        ny = nx
        n_level::Int = trunc(log(nx) / log(2.0) + 0.1)  # original c code uses natural log too
        h = xright / nx  # space step size defined as a global variable
        h2 = h^2 #space step size squared defined as a global variable
        if dt == 0
            dt = 0.1 * h2  # ∆t defined as a global variable
        end

        if gam == 0
            gam = M * h / (2 * sqrt(2) * atanh(0.9))
        end
        Cahn = gam^2  # ϵ^2 defined as a global variable

        #make directory
        if isdir(outdir)
            if !isempty(outdir)
                if overwrite == false
                    if check_dir == true
                        println("Warning: Directory is not empty. Results may be appended to existing files. Are you sure you want to continue? (Y/N)")
                        input = readline()
                        if input == "Y" || input == "y"
                            println("Appending to any existing files.")
                        else
                            println("End.")
                            break
                        end
                    end
                else
                    if check_dir == true
                        println("Warning: overwriting directory with new files. Are you sure you want to continue? (Y/N)")
                        input = readline()
                        if input == "Y" || input == "y"
                            rm(outdir, recursive=true)
                            mkdir(outdir)
                        else
                            println("End.")
                            break
                        end
                    end
                end
            end
        else
            mkpath(outdir)
        end

        #check if oc is in the range -1 to 1; if not, rescale to 0 to 1 and then shift with 2x - 1
        if minimum(oc) >= 0
            type = "psi"
            oc = (oc .- minimum(oc)) ./ (maximum(oc) - minimum(oc))
            oc = 2 .* oc .- 1    # psi0=(phi0+1)/2 
        else #if minimum(oc)>=-1 && maximum(oc)<=1
            type = "phi"
            # else
            # println("Warning: cannot interpret initial data. Data should be in range [-1,1] or only positive values.")
            # break
        end
        println(type)

        #initialization
        println("nx = $nx, ny = $ny, dt = $dt, epsilon = $gam, max_it = $max_it,max_it_CH= $max_it_CH, ns = $ns, n_level = $n_level")
        mu = zeros(Float64, nx, ny)
        nc = copy(oc)
        oc0 = copy(oc)
        if print_phi
            if type == "phi"
                open("$(outdir)/phi_$(nx)_$(max_it)_$(tol)_$(suffix).txt", "w", lock=false) do f
                    writedlm(f, oc, " ")
                end
            elseif type == "psi"
                open("$(outdir)/psi_$(nx)_$(max_it)_$(tol)_$(suffix).txt", "w", lock=false) do f
                    psi = (oc .+ 1) ./ 2
                    writedlm(f, psi, " ")
                end
            end
        end

        #run cahn hilliard
        for it in 1:max_it
            nc = cahn(oc, nc, mu, nx, ny, dt, max_it_CH, tol, c_relax, xright, xleft, Cahn, n_level, alpha, suffix=suffix)

            if print_mass
                print_mat("$(outdir)/ave_mass_$(nx)_$(max_it)_$(tol)_$(suffix).txt", calculate_mass(oc, h2, nx, ny))
            end
            if print_e
                print_mat("$(outdir)/discrete_norm_e_$(nx)_$(max_it)_$(tol)_$(suffix).txt", calculate_discrete_norm_energy(oc, oc0, h2, nx, ny, Cahn))
            end
            oc = copy(nc)
            if it % ns == 0
                if print_phi
                    println(it)

                    if type == "phi"
                        open("$(outdir)/phi_$(nx)_$(max_it)_$(tol)_$(suffix).txt", "a", lock=false) do f
                            writedlm(f, oc, " ")
                        end
                    elseif type == "psi"
                        open("$(outdir)/psi_$(nx)_$(max_it)_$(tol)_$(suffix).txt", "a", lock=false) do f
                            psi = (oc .+ 1) ./ 2
                            writedlm(f, psi, " ")
                        end
                    end
                end
            end
        end
        break
    end
end
