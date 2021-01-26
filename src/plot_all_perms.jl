using Plots
using Plots.PlotMeasures

include("initialize.jl")
include("atom_shuffle.jl")
include("sro_calc.jl")

file = "ab.yml" # Please note, that the file has to be in the same folder.
X = initialize(file)

all_perms = all_permutations(sequence(X))
all_sros = []

for i = 1:size(all_perms)[1]
    push!(all_sros, sro(X, all_perms[i]))
end

length(all_sros)
all_sros[1][2]
values = []
for i = 1:length(all_sros)
    push!(values, all_sros[i][1])
end
println(minimum(all_sros))
"""
ticks = ["AAAAAAAABBBBBBBB", "BBBBBBBBAAAAAAAA"]
x_axis = [i for i = 1:1:length(all_sros)]
length(values)
length(x_axis)
plt = plot(title = "All SRO values",
        xlabel = "Configuration",
        ylabel = "SRO value",
        xticks = (1:12869:12870, ticks),
        tickfontsize = 20,
        legend = false,
        palette = :bluesreds,
        bottom_margin = 10mm,
        left_margin = 10mm,
        right_margin = 35mm,
        size = (2000, 1300),
        guidefontsize = 25,
        titlefontsize = 35,
        )

# println(values)
plot!(plt, x_axis, values)
display(plt)
#savefig("all_sros.png")
"""
