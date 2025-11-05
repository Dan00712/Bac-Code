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
Z = range(start=0, stop=50e-6, length=200) |> x-> BigFloat.(x)
const κ = 2π *18e3

# function to minimize
function L(z, Δω, κ)
    f0 = α/ħ * Ec(0,0, z, Δω)^2
	δz = Δω - f0
	fz = im*k0 - im*zR/(z^2+zR^2) - z/(z^2+zR^2)
	gz = -2*z/Wc^2

	real(
		(1+f0/(δz+im*κ/2))*(fz+gz*f0/(δz-im*κ/2))
			)
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
         ylabel="z/μm",
         yaxis=:log,
         grid=true,
         minorgrid=true,
         formatter=:plain,
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

