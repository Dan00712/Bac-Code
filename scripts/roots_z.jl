using DrWatson
@quickactivate "SingleCavity"

include(scriptsdir("shared_code.jl"))
plotlyjs()
global_logger(ConsoleLogger(Info))

Ω = vcat(range(0, 400, length=200) .* 1e3*2π,  #.|> x-> 10^x,
         range(690e3, 700e3, length=50),
         range(740e3, 760e3, length=50)
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
         marker=:cross,
)

savefig(p, plotsdir("roots_z-$(now_nodots()).svg"))
gui(p)

if !Base.isinteractive()
    println("hit enter to close")
    readline()
end

