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

Ω = range(0, 400, length=200) .* 1e3*2π  #.|> x-> 10^x
#Ω = range(-10, ω0/1e6, length=75)
Z = range(start=-50e-6, stop=50e-6, length=200) |> x-> BigFloat.(x)
const κ = 2π *18e3

markers= [:rect, :star5, :xcross, :cross]
function L_carlos(z, Δ, κ)
	f0 = α/ħ * ħ*(Δ+ω0)/2/ϵ0/Vc * exp(-2*z^2/Wc^2)
    Δmod = Δ - f0

    (
        4 * z * (f0 + Δmod) * (2 * f0 * (z^2 + zR^2) + Wc^2 * Δmod)
        - 2 * f0 * Wc^2 * (-zR + k0 * (z^2 + zR^2)) * κ
        + Wc^2 * z * κ^2
       )/(Wc^2 * (z^2 + zR^2) * (4 * Δmod^2 + κ^2))
end
# function to minimize
function L(z, Δω, κ)
    L_carlos(z, Δω, κ)
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
         ylabel="z/μm",
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
         zmins.*1e6,
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

