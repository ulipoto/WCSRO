include("init.jl")
include("functions.jl")

function sro(Initializer, current_sequence)
    sequence = current_sequence
    X = Initializer

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
    for i = 1:X.total_atoms, j = i:X.total_atoms
        if X.distance_matrix[i, j] != 0
            X.pair_SRO[get_pair_key(X.pair_SRO, [sequence[i], sequence[j]])] += 1
            total_pairs += 1
        end
    end
    # Now, we calculate the true SRO parameter and store it in our dictionary
    # for pair SROs.
    for (key1, val1) in X.pair_SRO
        atoms = []
        for i = 1:2:length(key1)
            y = key1[i]*key1[i+1]
            push!(atoms, y)
        end
        X.pair_SRO[key1] = 1 - X.pair_SRO[key1]/((X.element_dict[atoms[1]]/X.total_atoms)*(X.element_dict[atoms[2]]/X.total_atoms)*X.total_atoms*12)
    end

    # Here, all triplet interactions are counted and stored in the triplet_SRO
    # dictionary.
    total_triplets = 0
    for i = 1:X.total_atoms, j = i:X.total_atoms, k = j:X.total_atoms
        if X.distance_matrix[i, j] != 0 && X.distance_matrix[i, k] != 0
            X.triplet_SRO[get_triplet_key(X.triplet_SRO, [sequence[i],
            sequence[j], sequence[k]])] += 1
            total_triplets += 1
        end
    end
    for (key1, val1) in X.triplet_SRO
        atoms = []
        for i = 1:2:length(key1)
            y = key1[i]*key1[i+1]
            push!(atoms, y)
        end
        X.triplet_SRO[key1] = 1 - (X.triplet_SRO[key1]/((X.element_dict[atoms[1]]/X.total_atoms)*
        (X.element_dict[atoms[2]]/X.total_atoms)*(X.element_dict[atoms[3]]/X.total_atoms)*X.total_atoms*12*11))
    end
    return X.pair_SRO, X.triplet_SRO
end
