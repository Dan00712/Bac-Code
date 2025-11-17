using DrWatson
@quickactivate "SingleCavity"

include(scriptsdir("shared_code.jl"))
global_logger(ConsoleLogger(Info))

Ω = vcat(range(0, 400, length=250) .* 1e3*2π,  #.|> x-> 10^x,
        ) |> sort
Z = range(start=-2-6, stop=log10(50)-6, length=1500) |> x-> 10 .^x #|> x-> BigFloat.(x)
Κ = 2π * [
          9e4,
          18e4,
          18e3
]

p = plot(;
#         title="κ=$((κ/2π/1000)) 2πkHz",
         xlabel="Δ/(2π kHz)",
         ylabel="z/μm",
         yaxis=:log,
         grid=true,
         minorgrid=true,
         formatter=:plain,
)

@info "generating data ∀ κ ∈ Κ"
for (κ, color, marker) in zip(Κ, ["orange", "green", "blue"] ,[:cross, :xcross, :star5])
    zmins = []
    ωs = []
    for ω in ProgressBar(Ω)
        newmins = get_zeroes(z->L(z, ω, κ), Z)
        append!(zmins, newmins)
        append!(ωs, [ω for _ in newmins])
    end
    
    scatter!(p,
             ωs/1e3,
             zmins.*1e6,
             label="κ=$(κ/2π /1e3) 2πkHz",
             marker=marker,
             color=color
    )
end

@info "save plot"

if Base.isinteractive()
    savefig(p, plotsdir("kappas-$(now_nodots()).png"))
    gui(p)
else
    savefig(p, plotsdir("kappas.tex"))
end

