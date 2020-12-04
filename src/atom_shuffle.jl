using Combinatorics
using Permutations
using Random

"""
shuffle_atoms(Initializer::Initializer, perms::Int)

This function creates 'perms' permutations of a given supercell.
"""
function shuffle_atoms(Initializer::Initializer, perms::Int)
    sequence = []
    for (key, val) in Initializer.element_dict
        for i = 1:val
            push!(sequence, key)
        end
    end
    unique_perms = Array{Any}(undef, perms, length(sequence)[1])
    unique_perms[1, :] = sequence
    i = 1
    while i < perms
        current_sequence = shuffle(sequence)
        v = 0
        for j = 1:i
            if current_sequence == unique_perms[j]
                v += 1
                break
            end
        end
        if v == 0
            i += 1
            unique_perms[i, :] = current_sequence
        end
    end
    return unique_perms
end

"""
sequence(X::Initializer)

A simple help function for creating all permutations.
"""
function sequence(X::Initializer)
    my_sequence = []
    for key in keys(X.element_dict)
        for i=1:X.element_dict[key]
            push!(my_sequence, key)
        end
    end
    return my_sequence
end

"""
all_permutations(x::T, prefix=T())

Generates all permutations of a given vector.
"""
function all_permutations(x::T, prefix=T()) where T
    if length(x) == 1
        return [[prefix; x]]
    else
        t = T[]
        for i in eachindex(x)
            if i > firstindex(x) && x[i] == x[i-1]
                continue
            end
            append!(t, all_permutations([x[begin:i-1];x[i+1:end]], [prefix; x[i]]))
        end
        return t
    end
end
