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


function isstable(x, y, z, Δω)
    isposdef(HL(x, y, z, Δω))
end


Δω = 100
X = (-zR):(2zR/10):zR .|> BigFloat
Y = (-zR):(2zR/10):zR .|> BigFloat
Z = ((-zR):(2zR/600):zR) ./ 1e15 .|> BigFloat


function stable_positions(Δω; 
        X=BigFloat.((-zR):(2zR/10):zR),
        Y=BigFloat.((-zR):(2zR/10):zR),
        Z=BigFloat.(((-zR):(2zR/200):zR) ./ 10)
    )
    Grid = Array{BigFloat}(undef, length(X), length(Y), length(Z))
    for ((i, x), (j, y), (k, z)) in
        (Iterators.product(enumerate(X), enumerate(Y), enumerate(Z)))
        Grid[i, j, k] = L(x, y, z, Δω)
    end
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
    
    for ((i, x), (j, y), (k, z)) in
        (Iterators.product(enumerate(X), enumerate(Y), enumerate(Z)))
        if isminima(Grid, i, j, k) && isstable(x,y,z,Δω)
            @debug "found minima at " i j k
            push!(Xm, x)
            push!(Ym, y)
            push!(Zm, z)
        end
    end

    Xm, Ym, Zm
end

lock = ReentrantLock()
Xm = BigFloat[]
Ym = BigFloat[]
Zm = BigFloat[]
Δωs = Float64[]

Threads.@threads for Δω in ProgressBar((-ω0:2ω0/40:ω0))
    Xm_, Ym_, Zm_ = stable_positions(Δω)

    @lock lock begin
        append!(Xm, Xm_)
        append!(Ym, Ym_)
        append!(Zm, Zm_)
        append!(Δωs, [Δω for _ in Xm])
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

scatter!(p,
         Xm*1e6, Ym*1e6, Zm*1e6
)

savefig(p, plotsdir("full_roots3d-omegas.now_nodots()).svg"))

gui(p)
