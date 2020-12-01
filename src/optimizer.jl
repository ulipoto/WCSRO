include("sro_calculator.jl")

function optimize(Initializer, permutations)
    """
    This function has the purpose to preserve some clarity in the WCSRO.jl file.
    Basically, here it is being iterated over all permutations that have been created
    in atom_shuffle.jl and the corresponding SRO parameter is being calculated.
    The two lowest SRO sums are returned. Note that the minima are treated separatedly.
    For now, there is no all-inclusive SRO sum due to lack of weighing factors etc.
    """
    X = Initializer
    perms = permutations
    best_pair_perm = 1000
    best_triplet_perm = 1000
    best_pair_config = perms[1]
    best_triplet_config = perms[1]
    for i = 1:size(perms)[1]
        pairs, triplets = sro(X, perms[1, :])
        current_pair_perm = 0
        for (idx, val) in pairs
            current_pair_perm += abs(val)
        end
        if current_pair_perm < best_pair_perm
            best_pair_config = perms[i, :]
        end
        current_triplet_perm = 0
        for (idx, val) in triplets
            current_triplet_perm += abs(val)
        end
        if current_triplet_perm < best_triplet_perm
            best_triplet_perm = current_triplet_perm
            best_triplet_config = perms[i, :]
        end
    end
    return best_pair_config, best_triplet_config
end
