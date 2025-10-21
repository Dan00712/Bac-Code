using DrWatson
@quickactivate "SingleCavity"
using LinearAlgebra
using Logging

global_logger(ConsoleLogger(Info))

using ProgressBars
using ForwardDiff
using Plots
plotlyjs()

using DelimitedFiles
using JLD2

using SingleCavity.Util
using SingleCavity.Constants
using SingleCavity.Laser

function L(x, y, z, Δω)
    -abs2(Et(x, y, z)) +
    abs2(αeq_c(x, y, z, Δω)*Ec(x, y, z, Δω)) +
    2*real(Et(x, y, z) * αeq_c(x, y, z, Δω) * Ec(x, y, z, Δω))
end

HL(x, y, z, Δω) = ForwardDiff.hessian(r->L(r..., Δω), [x, y, z])


function isminima(G, i, j, k)
    imax, jmax, kmax = size(G)

    if i == 1 || i == imax || j == 1 || j == jmax || k == 1 || k == kmax
        return false
    end

    if G[i, j, k] >= G[i+1, j, k] || G[i, j, k] >= G[i-1, j, k]
        return false
    end
    if G[i, j, k] >= G[i, j+1, k] || G[i, j, k] >= G[i, j-1, k]
        return false
    end

    if G[i, j, k] >= G[i, j, k+1] || G[i, j, k] >= G[i, j, k-1]
        return false
    end

    return true
end

function isstable(x, y, z, Δω)
    isposdef(HL(x, y, z, Δω))
end


Δω = 100
X = (-zR):(2zR/100):zR .|> BigFloat
Y = (-zR):(2zR/100):zR .|> BigFloat
Z = ((-zR):(2zR/200):zR) ./ 10 .|> BigFloat


@info "constructing Grid"
Grid = Array{BigFloat}(undef, length(X), length(Y), length(Z))
for ((i, x), (j, y), (k, z)) in
    ProgressBar(Iterators.product(enumerate(X), enumerate(Y), enumerate(Z)))
    Grid[i, j, k] = L(x, y, z, Δω)
end
@info "constructed Grid"
@debug "" G

if all(Grid .== 0)
    @error "all elements are zero"
    exit(1)
end
Xm = []
Ym = []
Zm = []
colors = []
stability = []

@info "checking for minima"
for ((i, x), (j, y), (k, z)) in
    ProgressBar(Iterators.product(enumerate(X), enumerate(Y), enumerate(Z)))
    if isminima(Grid, i, j, k)
        @debug "found minima at " i j k
        push!(Xm, x)
        push!(Ym, y)
        push!(Zm, z)
        push!(stability, if isstable(x, y, z, Δω)
            :stable
        else
            :unstable
        end)
    end
end

p = plot(;
    xlims = (X[1], X[end]) .* 1e6,
    xlabel = "x/μm",
    ylims = (Y[1], Y[end]) .* 1e6,
    ylabel = "y/μm",
    zlims = (Z[1], Z[end]) .* 1e6,
    zlabel = "z/μm",
)

group = map(stability) do state
    if state == :stable
        ""
    else
        "un"
    end * "stable equilibrium"
end

scatter!(
    p,
    Xm .* 1e6,
    Ym .* 1e6,
    Zm .* 1e6,
    group = group,
)

@info "writing outputs"
savefig(p, plotsdir("full_roots3d.N$(size(Grid))-$(now_nodots()).svg"))


@save datadir("sims", "roots_full", "sim.N$(size(Grid)).$(now_nodots()).jld2") Xm Ym Zm stability Grid

open(datadir("sims", "roots_full", "sim.N$(size(Grid))).$(now_nodots()).csv"), "w") do io
    writedlm(io, ["Xmins", "Ymins", "Zmins", "stability"], ",")
    writedlm(io, hcat(Xm, Ym, Zm, stability), ",")
end
gui(p)
