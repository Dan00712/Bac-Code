module Constants

export ħ, ϵ0, c,
        λ0, k0, ω0, Pt, Wt, Ax, Ay, E0, zR,
        ρ, ϵ, α, 
        Wc, Lc, Vc, κ, ζ

# using SI Units
const ħ = 6.626e-34       # Js
const ϵ0 = 8.854e-12      # F/m
const c = 3e+8            # m/s


# Tweezer

const λ0 = 1.55e-6        # m = 1.66μm
const k0 = 2π/λ0
const ω0 = c*k0

const Pt = 0.5            # W
const Wt = 1e-6           # m

const Ax = 1
const Ay = 1

const zR = 1/2 * k0 * Wt^2


const E0 = sqrt(4*Pt /π /ϵ0 /c /Wt^2 /Ax^2 /Ay)

# SiO₂
const ρ = 2200               # kg/m^3
const ϵ = 2.07
const α = (ϵ-1)ϵ/ρ


# Cavity
const Wc = 48e-6        # m
const Lc = 6.46e-3      # m
# const E0c =           # function of \omega_c
# keep ω_0 = const
const Vc = π*Wc^2*Lc/4

const κ = 2π * 1.06e6         # Hz
const ζ = 1

end # module
