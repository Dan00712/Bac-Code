module Newton

using ForwardDiff

export newton_1d, get_guesses

∂ = ForwardDiff.derivative
∇ = ForwardDiff.gradient


function get_guesses(f, range)
    xprev, prev = range[1], f(range[1])

    guesses = Vector{eltype(range)}()
    for x in range
        curr = f(x)
        if curr * prev < 0
            push!(guesses, (x + xprev)/2)
        end
        xprev = x
        prev = curr
    end
    guesses
end

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
