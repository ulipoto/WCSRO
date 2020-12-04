"""
sro(Initializer::Initializer, current_sequence::Array)

This function calculates an overall sum of both pair- and triplet-SRO parameters.
Note that there has to be weighing factors for both parts. This part, however,
is still under development. Therefore, sro() only provides a qualitative measure-
-ment for order.
"""
function sro(Initializer::Initializer, current_sequence::Array)
    sequence = current_sequence
    X = Initializer
    pair_SRO_sum = 0
    triplet_SRO_sum = 0

    # At first, the values in both dicts have to be cleared
    for (key, val) in X.pair_SRO
        X.pair_SRO[key] = 0
    end
    for (key, val) in X.triplet_SRO
        X.triplet_SRO[key] = 0
    end

    # Here, all pair interactions are counted and stored in the pair_SRO
    # dictionary.
    total_pairs = 0
    for i = 1:length(X.fcoords), j = i:length(X.fcoords)
        if X.distance_matrix[i, j] != 0 && X.distance_matrix[i, j] < 13
            X.pair_SRO[get_pair_key(X.pair_SRO, (sequence[i], sequence[j]))] += 1
            total_pairs += 1
        end
    end
    # Now, we calculate the true SRO parameter and store it in our dictionary
    # for pair SROs.
    for (key, val) in X.pair_SRO
        X.pair_SRO[key] = 1 - X.pair_SRO[key]/((X.element_dict[key[1]]/length(X.fcoords))*(X.element_dict[key[2]]/length(X.fcoords))*length(X.fcoords)*12)
        pair_SRO_sum += X.pair_SRO[key]
    end

    # Here, all triplet interactions are counted and stored in the triplet_SRO
    # dictionary.
    total_triplets = 0
    for i = 1:length(X.fcoords), j = i:length(X.fcoords), k = j:length(X.fcoords)
        if X.distance_matrix[i, j] != 0 && X.distance_matrix[i, k] != 0 && X.distance_matrix[i, j] < 13 && X.distance_matrix[i, k] < 13
            X.triplet_SRO[get_triplet_key(X.triplet_SRO, (sequence[i],
            sequence[j], sequence[k]))] += 1
            total_triplets += 1
        end
    end
    for (key, val) in X.triplet_SRO
        X.triplet_SRO[key] = 1 - (X.triplet_SRO[key]/((X.element_dict[key[1]]/length(X.fcoords))*
        (X.element_dict[key[2]]/length(X.fcoords))*(X.element_dict[key[3]]/length(X.fcoords))*length(X.fcoords)*12*11))
        triplet_SRO_sum += X.triplet_SRO[key]
    end
    sro = pair_SRO_sum + triplet_SRO_sum
    return sro
end

"""
optimize(Initializer, permutations)

This function has the purpose to preserve some clarity in the WCSRO.jl file.
However, this is not used at the moment to have some freedom with creating plots etc.
Basically, here it is being iterated over all permutations that have been created
in atom_shuffle.jl and the corresponding SRO parameter is being calculated.
The two lowest SRO sums are returned. Note that the minima are treated separatedly.
For now, there is no all-inclusive SRO sum due to lack of weighing factors etc.
"""
function optimize(Initializer, permutations)
    X = Initializer
    perms = permutations
    best_pair_perm = 1000
    best_triplet_perm = 1000
    best_pair_config = perms[1]
    best_triplet_config = perms[1]
    for i = 1:size(perms)[1]
        pairs, triplets = sro(X, perms[1, :])
        println(pairs)
        current_pair_perm = 0
        for (idx, val) in pairs
            println(abs(val))
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
