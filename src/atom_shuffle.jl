using Combinatorics
using Permutations
using Random

include("init.jl")
include("structs.jl")
function shuffle_atoms(sequence, perms)
    sequence = sequence
    unique_perms = []
    push!(unique_perms, sequence)

    while length(unique_perms) < perms
        current_seq = shuffle(sequence)
        v = 0
        for n = 1:length(unique_perms)
            if current_seq ==  unique_perms[n]
                v += 1
                break
            end
        end
        if v == 0
            push!(unique_perms, current_seq)
        end
    end
    return unique_perms
end
