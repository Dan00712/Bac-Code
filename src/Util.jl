module Util

using Dates

export now_nodots, isminima

function now_nodots()
    now() |> string |> x->split(x, ".")[1] |> x -> replace(x, ":" => "-")
end

function isminima(G, i, j, k)
    imax, jmax, kmax = size(G)

    if i == 1 || i == imax || j == 1 || j == jmax || k == 1 || k == kmax
        return false
    end

    if G[i, j, k] >= G[i+1, j, k] || G[i, j, k] >= G[i-1, j, k]
        return false
    end
    if G[i, j, k] >= G[i, j+1, k] || G[i, j, k] >= G[i, j-1, k]
        return false
    end

    if G[i, j, k] >= G[i, j, k+1] || G[i, j, k] >= G[i, j, k-1]
        return false
    end

    return true
end

end #module
