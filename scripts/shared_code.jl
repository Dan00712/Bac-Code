using LinearAlgebra
using Logging

using ProgressBars

using Plots
if isinteractive()
    plotlyjs()
else
    gr()
end

using DelimitedFiles
using JLD2

using SingleCavity.Util
using SingleCavity.Constants
using SingleCavity.Laser

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

function isstable(z, Δω, κ)
	S = HHs(0,0, z, Δω, κ)
	N = Int(size(S, 1) // 2)
	J = let
		A = I(N) .|> ComplexF64
		A[N, N] = 1/im /ħ
		[zeros(N,N)  A ; -A zeros(N,N)]
	end
	K = let
		B = zeros(N, N) .|> ComplexF64
		B[N,N] = κ/2
		[
			B 			zeros(N,N);
			zeros(N,N) 	B
		]
	end

	A = J*S - K
	V = eigen(A).values
	ReV = real.(V)

    tol = 1e-9
	#@show ma, mi
	foo(x) = isapprox(x, 0; atol=tol) || x < 0
	foo.(ReV) |> all
end
