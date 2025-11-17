using LinearAlgebra
using Logging

using ProgressBars

using Plots
if isinteractive()
    plotlyjs()
else
    pgfplotsx()
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

function issemiposdef(A; tol = 1e-10)
    H = Hermitian(A)  # ensures Hermitian view, symmetrizes if needed
    λmin = minimum(eigvals(H))
    return λmin >= -tol
end

isstable(z, Δω, κ) = isposdef(HHs(0, 0, z, Δω, κ))
issemistable(z, Δω, κ) = issemiposdef(HHs(0, 0, z, Δω, κ))


