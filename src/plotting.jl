using Plots
using Plots.PlotMeasures
using Random
using LaTeXStrings
using StatsPlots
include("shuffle_config.jl")
include("initialize.jl")
include("sro_calc.jl")


"""
plot_configs(my_sros, filename="my_config_plot.png")


"""
function plot_configs(my_sros, filename="my_config_plot.png")
    #ticks = ["aaaaaaaabbbbbbbb", "bbbbbbbbaaaaaaaa"]
    x_axis = [i for i = 1:1:length(my_sros)]
    my_plt = plot(title = "SRO values",
            xlabel = "Configuration",
            ylabel = "Objective function",
            xticks = (1:length(my_sros)-1:length(my_sros)),#, ticks),
            tickfontsize = 20,
            legend = false,
            palette = :rainbow,
            bottom_margin = 10mm,
            left_margin = 10mm,
            right_margin = 35mm,
            size = (2000, 1300),
            guidefontsize = 25,
            titlefontsize = 35,
            )

    plot!(my_plt, x_axis, my_sros)
    display(my_plt)
end

"""
xyzplot(Backpack, best_sequence)

Just a little script for displaying the minimum structure. If you do not like my
chosen color palette, feel free to change it ;-)
"""
function xyzplot(X::Backpack, seq, cam_angle::Float64)
    plt = scatter(title = "Supercell", xlabel = "X in Å", ylabel = "Y in Å", gridlinewidth = 1, zlabel = "Z in Å", camera = (cam_angle, 15))
    my_colors = palette(:Accent_8) #Set3_12
    k = 0
    for (idx, val) in X.element_dict
        xyz = zeros(length(seq), 3)
        j = 0
        k += 1
        for i = 1:length(seq)
            if seq[i] == val
                j += 1
                xyz[i, :] =  X.supercell.lattice*X.supercell.fcoords[i, :]
            end
        end
        scatter!(plt, xyz[:, 1], xyz[:, 2], xyz[:, 3], markersize = 11, markercolor = my_colors[k],  label=idx)
    end
    return(plt)
end


"""
xyz_animation(X::Backpack, seq, gifname="supercell.gif", fps=15)

Animates a ittle gif of the supercell for better viewability
"""
function xyz_animation(X::Backpack, seq, gifname="supercell.gif", fps=15)
    n = 100
    t1 = [i*90/n for i=1:n]
    t = vcat(t1, reverse(t1))
    anim = @animate for i=1:length(t)
        xyzplot(X, seq, t[i])
    end
    return gif(anim, gifname, fps = 15)
end
