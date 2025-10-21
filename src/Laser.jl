module Laser

using ..Constants

export Ec, Et, αeq_c

function E0c(Δω)
    # Δω = ωc - ω0
    ωc = Δω + ω0

    sqrt(ħ*ωc/2/ϵ0/Vc)
end


function Ec(x, y, z, Δω)
    kc = let
        ωc = Δω + ω0
        ωc/c
    end
    E0c(Δω) * exp(-(x^2+z^2)/Wc^2) * ζ * cos(kc)
end

function ϕt(x, y, z)
    -atan(z/zR) + k0*z/2 * (x^2+y^2)/(z^2+zR^2)
end

function W(z)
    Wt * sqrt(1+(z/zR)^2)
end

expi(z) = exp(im*z)

function Et(x, y, z)
    E0/2 * expi(k0*z+ϕt(x, y, z)) * Wt/W(z) * exp(-(y/Ay/W(z))^2) * exp(-(x/Ax/W(z)^2)^2)
end

function δc(x, y, z, Δω)
    Δω - α/ħ * (Ec(x, y, z, Δω))^2
end

function αeq_c(x, y, z, Δω)
    α/ħ * Et(x, y, z) * Ec(x, y, z, Δω)/(δc(x, y, z, Δω)-im*κ/2)
end

end # module
