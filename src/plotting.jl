using Plots
using Random

function xyzplot(Initializer, best_sequence)
    """
    Just a little script for displaying the minimum structure. If you hate my
    chosen color palette, feel free to change it ;-)
    """
    X = Initializer
    sequence = best_sequence
    plt = scatter()
    my_colors = palette(:Set3_12)
    k = 0
    for (idx, val) in X.element_dict
        j = 0
        k += 1
        xyz = zeros(val, 3)
        for i = 1:length(sequence)
            if sequence[i] == idx
                j += 1
                xyz[j, 1] = X.xyz[i, 1]
                xyz[j, 2] = X.xyz[i, 2]
                xyz[j, 3] = X.xyz[i, 3]
            end
        end
        scatter!(plt, xyz[:, 1], xyz[:, 2], xyz[:, 3], markersize = 11, color = my_colors[k],  label=idx)
    end
    display(plt)
end
