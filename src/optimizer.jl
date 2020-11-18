include("sro_calculator.jl")

function optimize(Initializer, permutations)

    """
    This function has the purpose to preserve some clarity in the WCSRO.jl file.
    Basically, here one iterates through all permutations that have been created
    in atom_shuffle.jl and calculate the corresponding SRO parameter. The two
    lowest SRO sums are returned. Note that the minima are treated separatedly.
    For now, there is no all-inclusive SRO sum due to lack of weighing factors etc.
    """
    X = Initializer
    perms = permutations
    best_pair_perm = 1000
    best_triplet_perm = 1000
    best_pair_config = my_permutations[1]
    best_triplet_config = my_permutations[1]
    for i = 1:length(perms)
        pairs, triplets = sro(X, my_permutations[i])
        current_pair_perm = 0
        for (idx, val) in pairs
            current_pair_perm += abs(val)
        end
        #push!(pair_perms, current_pair_perm)
        if current_pair_perm < best_pair_perm
            best_pair_config = my_permutations[i]
        end
        current_triplet_perm = 0
        for (idx, val) in triplets
            current_triplet_perm += abs(val)
        end
        #push!(triplet_perms, current_triplet_perm)
        if current_triplet_perm < best_triplet_perm
            best_triplet_perm = current_triplet_perm
            best_triplet_config = my_permutations[i]
        end
    end
    return best_pair_config, best_triplet_config
end
