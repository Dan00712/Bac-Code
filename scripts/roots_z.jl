using DrWatson
@quickactivate "SingleCavity"
using LinearAlgebra
using Logging


using ProgressBars

using ForwardDiff
∂ = ForwardDiff.derivative

using Plots
plotlyjs()

using DelimitedFiles
using JLD2

using SingleCavity.Util
using SingleCavity.Constants
using SingleCavity.Laser
using SingleCavity.Newton


# We can optimize the method used to a 1d root finding problem
# since it can be analytically shown that x=y=0 is a requirement we only need to do this for a single function
# but for checking the full hamiltonian is needed, so this will only be computed on a as needed basis

# to determine guesss for roots we will use the method of creating samples and checking for sign changes
# this will miss any zero positions where the root only touches zero

global_logger(ConsoleLogger(Info))

function δc(x,y,z, Δω)
    Δω - α/ħ * Ec(x,y,z, Δω)^2
end

function γ(x,y,z, Δω)
    α/ħ * Ec(x,y,z, Δω)/(δc(x,y,z, Δω)^2 + κ^2/4)
end

function Hs(x, y, z, Δω)
    -abs2(Et(x, y, z)) +
    abs2(αeq_c(x, y, z, Δω)*Ec(x, y, z, Δω)) +
    2*real(Et(x, y, z) * αeq_c(x, y, z, Δω) * Ec(x, y, z, Δω))
end

# Hessian of the Hamiltonian
function HHs(x,y,z, Δω)
    ForwardDiff.hessian(r->Hs(r..., Δω), [x, y, z])
end

# function to minimize
function L(z, Δω)
    -(z * 
        (1+ γ(0,0,z, Δω)*δc(0,0,z, Δω)) * 
        (1/(zR^2 + z^2) + 2/Wc^2 * γ(0,0,z, Δω)*δc(0,0,z, Δω)) +
    γ(0,0,z, Δω) * κ/2 * 
    (2*z/Wc^2 + zR/(z^2+zR^2) - k0)) * abs2(Et(0,0, z))
end

# Δω = 0:2ω0
Ωs = -ω0:2ω0/100:ω0
Zs = (-zR:2zR/200:zR) ./1e2
@debug "" Zs length(Zs)

function get_roots(Δω)
    f(z) = L(z, Δω)

    guesses = get_guesses(f, Zs)

    zmin = guesses

    zmin
end

ωs    = BigFloat[]
zmins = BigFloat[]


@info "Δω = 0" get_minimas(0)

p = plot()

for  ω in Ωs
    zmins_ = get_minimas(ω)
    @info "" zmins_

    append!(zmins, zmins_)
    append!(zmins, [ω for _ in zmins_])
end
scatter!(p, ωs, zmins)
#plot!(p, Zs, L.(Zs,[0]))
#plot!(p, Zs, (z->∂(x-> Hs(0,0,x, 0), z)).(Zs))

gui(p)


