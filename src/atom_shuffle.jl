using Combinatorics
using Permutations
using Random
include("init.jl")

function shuffle_atoms(perms)
    global unique_perms = []
    push!(unique_perms, sequence)
    total_permutations = factorial(big(a_counter + b_counter + c_counter))/(factorial(a_counter)*factorial(b_counter)*factorial(c_counter))

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
end
