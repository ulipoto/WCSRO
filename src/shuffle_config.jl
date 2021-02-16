using Random
include("structure.jl")
include("initialize.jl")

"""
shuffle_atoms(Backpack::Backpack, perms::Int)

This function creates 'perms' permutations of a given supercell.
"""
function shuffle_atoms(Backpack::Backpack, perms::Int)
    X = Backpack
    seq = [X.element_dict[key] for key in X.supercell.species]
    unique_perms = Array{Vector{Tuple{Int64,Int64}}}(undef, perms)
    unique_perms[1] = seq
    i = 1
    while i < perms
        current_sequence = sublattice_shuffle(seq, length(X.supercell.composition))
        v = 0
        for j = 1:i
            if current_sequence == unique_perms[j]
                v += 1
                break
            end
        end
        if v == 0
            i += 1
            unique_perms[i] = current_sequence
        end
    end
    return unique_perms
end
function sublattice_shuffle(seq, sublattices::Int)
    new_seq = Array{Any}(undef, length(seq))
    for i = 1:sublattices
        indices = findall(x -> x[1] == i, seq)
        new_indices = shuffle(indices)
        for j = 1:length(indices)
            new_seq[indices[j]] = seq[new_indices[j]]
        end
    end
    return new_seq
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
