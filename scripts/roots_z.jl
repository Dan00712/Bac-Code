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

global_logger(ConsoleLogger(Info))

Ω = range(-10, 1e10, length=50)
Z = range(start=-5e2zR, stop=5e2zR, length=100) |> x-> BigFloat.(x)
const κ = 2π *1.06e6/1e1

# function to minimize
function L(z, Δω, κ)
    fz = im * (k0 - zR/(z^2+zR^2)) - z/(z^2+zR^2)
    gz = -2*z/Wc^2
    γz = γ(0,0,z, Δω, κ)
    δz = δc(0,0,z,Δω)

    (1+γz*δz)*(real(fz)+ gz*γz*δz) + γz * κ/2 * (imag(fz) + gz*γz*κ/2)
end

function get_zeros(f)
    zmins = []

    for z in Z
        try
            zmin = find_zero(f, z)
            if zmin ∉ zmins
                push!(zmins, zmin)
            end
        catch e
            if !isa(e, Roots.ConvergenceFailed)
                rethrow(e)
            end
        end
    end

    zmins
end

p = plot(;
         xlabel="Δ/2π",
         #xaxis=:log,
         ylabel="z/zR",
)
zmins = []
ωs = []
stabilities = []
for ω in ProgressBar(Ω)
    newmins = get_zeros(z->L(z, ω, κ))
    append!(zmins, newmins)
    append!(ωs, [ω for _ in newmins])
end

for (ω, z) in zip(ωs, Z)
    push!(stabilities,
          isposdef(HHs(0,0,z, ω, κ)) ? :stable : :unstable)
end

scatter!(p,
         ωs,
         zmins./zR,
         label="κ = $(round(κ/2π)) 2πHz",
         group=stabilities,
)
savefig(p, plotsdir("full_z-$(now_nodots()).svg"))
gui(p)

if !Base.isinteractive()
    println("hit enter to close")
    readline()
end

