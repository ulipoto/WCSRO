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

    for i = 1:length(sequence), j = i+1:length(sequence)
        if X.shell_matrix[i, j] != 0
            X.pair_SRO[X.element_dict[sequence[i]][1], X.element_dict[sequence[j]][1], X.shell_matrix[i, j]] +=
            X.pair_weight_matrix[X.element_dict[sequence[i]][1], X.element_dict[sequence[j]][1], X.shell_matrix[i, j]]
        end
    end

    # Here, all triplet interactions are counted and stored in the triplet_SRO
    # dictionary.
    for i = 1:size(X.shell_matrix)[1], j = i:size(X.shell_matrix)[1], k = j:size(X.shell_matrix)[1]
        if X.shell_matrix[i, j] != 0
            X.triplet_SRO[X.element_dict[sequence[i]][1], X.element_dict[sequence[j]][1], X.element_dict[sequence[k]][1], (X.shell_matrix[i, j]+X.shell_matrix[i, k])] +=
            X.triplet_weight_matrix[X.element_dict[sequence[i]][1], X.element_dict[sequence[j]][1], X.element_dict[sequence[k]][1], (X.shell_matrix[i, j]+X.shell_matrix[i, k])]
        end
    end

    return X.pair_SRO, X.triplet_SRO
end
