using Plots

include("init.jl")
include("sro_calculator.jl")

function xyzplot(best_sequence)

    sequence = best_sequence
    x1 = []
    y1 = []
    z1 = []
    x2 = []
    y2 = []
    z2 = []
    x3 = []
    y3 = []
    z3 = []

    for n = 1:length(sequence)
        if sequence[n] == atom_a
            push!(x1, xyz[n, 1])
            push!(y1, xyz[n, 2])
            push!(z1, xyz[n, 3])
        elseif sequence[n] == atom_b
            push!(x2, xyz[n, 1])
            push!(y2, xyz[n, 2])
            push!(z2, xyz[n, 3])
        else
            push!(x3, xyz[n, 1])
            push!(y3, xyz[n, 2])
            push!(z3, xyz[n, 3])
        end
    end
    plt = Plots.scatter3d(x1, y1, z1, markersize = 11, color= [:blue], label=[atom_a], title = "Current structure")
    Plots.scatter3d!(x2, y2, z2, markersize = 11, color= [:red], label=[atom_b])
    display(plt)
    if c_counter != 0
    Plots.scatter3d!(x3, y3, z3, markersize = 11, color= [:green], label=[atom_c])
    end
end
