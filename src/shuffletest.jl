using IterTools
using LinearAlgebra
using Combinatorics
using Plots

include("structs.jl")
include("initialize.jl")
include("atom_shuffle.jl")
include("sro_calculator.jl")
include("plotting.jl")

file = "POSCAR.mp-1184260_Fe3Cu"

function unique_permutations(x::T, prefix=T()) where T
    if length(x) == 1
        return [[prefix; x]]
    else
        t = T[]
        for i in eachindex(x)
            if i > firstindex(x) && x[i] == x[i-1]
                continue
            end
            append!(t, unique_permutations([x[begin:i-1];x[i+1:end]], [prefix; x[i]]))
        end
        return t
    end
end

#X = initialize(file)
function generate_sequence(X::Initializer)
    sequence = []
    for key in keys(X.element_dict)
        for i=1:X.element_dict[key]
            push!(sequence, key)
        end
    end
    return sequence
end
println(generate_sequence(X))


function unique_permutations(x::T, prefix=T()) where T
    if length(x) == 1
        return [[prefix; x]]
    else
        t = T[]
        for i in eachindex(x)
            if i > firstindex(x) && x[i] == x[i-1]
                continue
            end
            append!(t, unique_permutations([x[begin:i-1];x[i+1:end]], [prefix; x[i]]))
        end
        return t
    end
end
z = ["A", "B", "C", "D"]
println(unique_permutations(z))
println(length(unique_permutations(z)))
"""
sro_list = []
for i = 1:length(all_perms)
    push!(sro_list, sro(X, all_perms[i]))
end

x_axis = []
for i = 1:size(all_perms)[1]
    current_sequence = ""
    for j = 1:length(all_perms[1])
        current_sequence = current_sequence * all_perms[i][j]
    end
    push!(x_axis, current_sequence)
end

k = plot(x_axis, sro_list,
    #seriestype = :scatter,
    size = (1500, 750),
    gridlinewidth = 2,
    label = "Supercell Fe12Cu4",
    left_margin = 22Plots.mm,
    right_margin = 50Plots.mm,
    bottom_margin = 10Plots.mm,
    title = "SRO parameters",
    xlabel = "Configurational space",
    ylabel = "SRO value",
    color = :blue4,
    xtickfont=font(14),
    ytickfont=font(14),
    guidefont=font(18),
    legendfont=font(14),
    titlefont=font(28),
    xticks = 3,
    savefig("plot.png")
    )
display(k)
"""
