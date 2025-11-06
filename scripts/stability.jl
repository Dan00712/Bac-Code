using DrWatson
@quickactivate "SingleCavity"

include(scriptsdir("shared_code.jl"))
pgfplotsx()

global_logger(ConsoleLogger(Info))

Ω = vcat(range(0, 400, length=100) .* 1e3*2π,  #.|> x-> 10^x,
         range(690e3, 760e3, length=80)
        ) |> sort
Z = range(start=-2-6, stop=log10(50)-6, length=1500) |> x-> 10 .^x #|> x-> BigFloat.(x)
const κ = 2π *18e3

p = plot(;
         title="κ=$((κ/2π/1000)) 2πkHz",
         xlabel="Δ/(2π kHz)",
         ylabel="z/μm",
         yaxis=:log,
         grid=true,
         minorgrid=true,
         formatter=:plain,
)

@info "checking equilibrium positions ∀ ω ∈ Ω"
zmins = []
ωs = []
for ω in ProgressBar(Ω)
    newmins = get_zeroes(z->L(z, ω, κ), Z)
    append!(zmins, newmins)
    append!(ωs, [ω for _ in newmins])
end


@info "check stability of equilibrium positions"
stabilities = [ isstable(zmin, ω, κ) ? :stable : :unstable
               for (ω, zmin) in zip(ωs, zmins)
              ]

@info "save generated data"
date = now_nodots()
M = (Δω=ωs, zmin=zmins, stability=stabilities)
mkpath(datadir("sims", "stability"))
@save datadir("sims", "stability", "s-$date.jld2") M

@info "do plot"
for (g, marker, color) in zip([:stable, :unstable], [:cross, :xcross], [:orange, :blue])
    mask = stabilities .== g
    scatter!(p,
            ωs[mask]/1e3,
            zmins[mask].*1e6,
            marker=marker,
            color=color,
            label=string(g),
    )
end

@info "save plot"
savefig(p, plotsdir("stabilities-$date.tex"))
savefig(p, plotsdir("stabilities.tex"))

if Base.isinteractive()
    gui(p)
end

