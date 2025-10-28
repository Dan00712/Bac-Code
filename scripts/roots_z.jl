using DrWatson
@quickactivate "SingleCavity"
using LinearAlgebra
using Logging

using ProgressBars

using Roots

using Plots
plotlyjs()

using DelimitedFiles
using JLD2

using SingleCavity.Util
using SingleCavity.Constants
using SingleCavity.Laser
using SingleCavity.Newton


global_logger(ConsoleLogger(Info))

# function to minimize
function L(z, Δω, κ)
    ∂zHs(0,0, z, Δω, κ)
end


Ω = range(start=0, stop=10e6, length=100)
Z = range(start=-5zR, stop=5zR, length=200)
const κ = 2π * 1.06e6 /1e1


function get_zeros(f)
    zmin = unique(find_zero(f, z) for z in Z)

    zmin
end


zmins = []
ωs = []
for ω in Ω
    newmins = get_zeros(z->L(z, ω, κ))
    append!(zmins, newmins)
    append!(ωs, [ω for _ in newmins])
end

p = plot(;
         xlabel="Δ/2π",
         ylabel="z/zR"
)

scatter!(p, ωs, zmins./zR)

savefig(p, plotsdir("full_z-$(now_nodots()).svg"))
gui(p)

if !Base.isinteractive()
    println("hit enter to close")
    readline()
end
