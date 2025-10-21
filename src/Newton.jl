module Newton

using ForwardDiff

∂ = ForwardDiff.derivative
∇ = ForwardDiff.gradient

function newton_1d(init, f::Function; max_iters = 100, atol = 1e-10)
    x = init
    xprev = init
    df(x) = ∂(f, x)

    for _ = 1:max_iters
        xprev = x
        x = x - f(x)/df(x)

        if abs(x - xprev) < atol
            break
        end
    end

    convergent = abs(x-xprev) < atol

    x, convergent
end

end # module
