using DrWatson
@quickactivate "SingleCavity"
include(scriptsdir("shared_code.jl"))
using ForwardDiff


function H(r, p, Δ, κ)
    x,y,z, ar, ai = r...
    px,py,pz, ars, ais = p...
    ωc = Δ + ω0

    norm(p)^2/2/m + ħ*ωc * (ar+im*ai) * (ars + im*ais) - α/2 * abs2(Et(x,y,z) + a * Ec(x,y,z, Δ))
end
∇H(r,p, Δ, κ) = ForwardDiff.gradient(f -> H(f[1:end/2], f[end/2:end], Δ, κ), [r..., p...])

function step(r, p, Δ, κ; dt)
    dH = ∇H(r, p, Δ, κ) * dt

    # rescaling for the Poisson bracket thing
    dH[end/2] *= i/ħ
    dH[end/2-1] *= i/ħ
    dH[end] *= i/ħ
    dH[end-1] *= i/ħ

    r + dH[1:end/2], p + dH[end/2:end]
end
