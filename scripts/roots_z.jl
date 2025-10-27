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
    γδ = γ(0,0,z, Δω) * δc(0,0,z, Δω)

    z * (1+γδ)*(1/zR^2 + 2/Wc^2 * γδ) +
        γ(0,0, z, Δω) * κ/2 * (2*z/Wc^2 + 1/zR - k0)
end

# Δω = 0:2ω0
Ωs = 10 .^(-2:12/200:10)
Zs = (-zR:2zR/1000:zR) |> x->BigFloat.(x)
#Zs = range(-zR/5, zR/5; length=200)
@debug "" Zs length(Zs)

function get_roots(Δω; atol=1e-10)
    f(z) = L(z, Δω)

    guesses = get_guesses(f, Zs)

    zmin = BigFloat[]
    for guess in guesses
        x, convergent = newton_1d(guess, f, atol=atol)
        if convergent
            push!(
                  zmin,
                  x
                 )
        end
    end

    zmin
end

ωs    = Float64[]; #BigFloat[]
zmins = Float64[]; #BigFloat[]

#@info "Δω = 0" get_roots(0)


for  ω in ProgressBar(Ωs)
    zmins_ = get_roots(ω, atol=1e-30)
    #@info "" zmins_

    append!(zmins, zmins_)
    append!(ωs, [ω for _ in zmins_])
end
p = plot(;
         #ylims=(-1, 1)./1e32,
         xlabel="Δω/[s^-1]",
         ylabel="z/[m]",
         xaxis=:log,
        )

scatter!(p, ωs, zmins)
#plot!(p, Zs, L.(Zs,[0]))
#plot!(p, Zs, (z->∂(x-> Hs(0,0,x, 0), z)).(Zs))
savefig(p, plotsdir("full_z.N$(size(Zs))-$(now_nodots()).svg"))

gui(p)


