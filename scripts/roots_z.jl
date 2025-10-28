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
    γδ = γ(0,0,z, Δω, κ) * δc(0,0,z, Δω)

    z * (1+γδ)*(1/zR^2 + 2/Wc^2 * γδ) +
        γ(0,0, z, Δω, κ) * κ/2 * (2*z/Wc^2 + 1/zR - k0)
end


Ω = range(start=0, stop=10e6, length=200)
Z = range(start=-5zR, stop=5zR, length=800)
const κ = 2π * 1.06e6 /1e2


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
)

scatter!(p, ωs, zmins)

savefig(p, plotsdir("full_z-$(now_nodots()).svg"))
gui(p)

if !Base.isinteractive()
    println("hit enter to close")
    readline()
end
