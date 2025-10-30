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
    fz = im * (k0 - zR/(z^2+zR^2)) - z/(z^2+zR^2)
    gz = -2*z/Wc^2
    γz = γ(0,0,z, Δω, κ)
    δz = δc(0,0,z,Δω)

    (1+γz*δz)*(real(fz)+ gz*γz*δz) + γz * κ/2 * (imag(fz) + gz*γz*κ/2)
end

Ω = range(0, 1e6, length=75) #(range(start=-5, stop=8, length=75) .|> x->10^x)
Z = range(start=-25zR, stop=25zR, length=100) |> x-> BigFloat.(x)
const κ = 2π *1.06e6/10


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
for ω in ProgressBar(Ω)
    newmins = get_zeros(z->L(z, ω, κ))
    append!(zmins, newmins)
    append!(ωs, [ω for _ in newmins])
end

scatter!(p,
         ωs,
         zmins./zR,
         label="κ = $(round(κ/2π)) 2πHz"
)
savefig(p, plotsdir("full_z-$(now_nodots()).svg"))
gui(p)

if !Base.isinteractive()
    println("hit enter to close")
    readline()
end

