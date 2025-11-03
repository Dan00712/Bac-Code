using DrWatson
@quickactivate "SingleCavity"
using LinearAlgebra
using Logging

using ProgressBars

using Roots

using Plots
plotlyjs()
#gr()

using DelimitedFiles
using JLD2

using SingleCavity.Util
using SingleCavity.Constants
using SingleCavity.Laser

global_logger(ConsoleLogger(Info))

Ω = range(0, 400e3, length=200) #.|> x-> 10^x
#Ω = range(-10, ω0/1e6, length=75)
Z = range(start=-2e2zR, stop=2e2zR, length=75) |> x-> BigFloat.(x)
const κ = 2π *1.06e6/5e2
#=Κ = 2π*1.06 *[
                1e6,
                1e5,
                5e5,
                1e4,
       ]=#
markers= [:rect, :star5, :xcross, :cross]
# function to minimize
function L(z, Δω, κ)
    fz = im * (k0 - zR/(z^2+zR^2)) - z/(z^2+zR^2)
    gz = -2*z/Wc^2
    γz = γ(0, 0, z, Δω, κ)
    δz = δc(0, 0, z, Δω)

    (1+γz*δz)*(real(fz)+ gz*γz*δz) + γz * κ/2 * (imag(fz) + gz*γz*κ/2)
end

function get_zeros(f)
    	guesses = let
		guesses = []
		zprev = Z[1]
		fprev = f(Z[1])

		for z in Z[2:end]
			fz = f(z)

			if fz * fprev < 0
                push!(guesses, (zprev, z))
			end
			zprev = z
			fprev = fz
		end
		guesses
	end

	zmins = []
	for guess in guesses
		push!(zmins, find_zero(f, guess))
	end

    zmins
end

p = plot(;
         xlabel="Δ/2π",
#         xaxis=:log,
         ylabel="z/zR",
#         legend=:bottom,
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
#         label="κ = $((κ/2π)) 2πHz",
#         group=stabilities,
marker=:cross,
)

savefig(p, plotsdir("full_z-$(now_nodots()).svg"))
gui(p)

if !Base.isinteractive()
    println("hit enter to close")
    readline()
end

