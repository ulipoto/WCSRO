using Random
using Combinatorics
include("atom_shuffle.jl")
include("initialize.jl")

"""
sro(Backpack::Backpack, current_sequence::Array)

This function calculates an overall sum of both pair- and triplet-SRO parameters.
A weight is integrated with 1/i, where i is the coordination shell number.
"""
function sro(Backpack::Backpack, sequence)
    X = Backpack
    seq = sequence
    pair_SRO_sum = zeros(size(X.pair_weight_matrix)[3])
    triplet_SRO_sum = zeros(size(X.triplet_weight_matrix)[4])
    P = copy(X.pair_SRO)
    T = copy(X.triplet_SRO)
    # Pair part
    for i = 1:length(X.all_pairs)
        X.pair_SRO[seq[X.all_pairs[i][1]], seq[X.all_pairs[i][2]], X.all_pairs[i][3]] +=
            X.pair_weight_matrix[seq[X.all_pairs[i][1]], seq[X.all_pairs[i][2]], X.all_pairs[i][3]]
    end
    for i = 1:length(pair_SRO_sum)
        for j = 1:length(X.element_dict), k = j:length(X.element_dict)
            if j != k
                pair_SRO_sum[i] += 1 - (X.pair_SRO[j, k, i] + X.pair_SRO[k, j, i])
            end
        end
    end
    X.pair_SRO = P
    total_pair_sum = 0
    for i = length(pair_SRO_sum)
        total_pair_sum += 1/i*(pair_SRO_sum[i])
    end

    # Triplet part
    for i = 1:length(X.all_triplets)
        idx = [seq[X.all_triplets[i][1]], seq[X.all_triplets[i][2]], seq[X.all_triplets[i][3]]]
        shell = X.all_triplets[i][4]
        sort!(idx)
        X.triplet_SRO[idx[1], idx[2], idx[3], shell] += X.triplet_weight_matrix[idx[1], idx[2], idx[3], shell]
    end
    for i = 1:length(triplet_SRO_sum)
        for j = 1:length(X.element_dict), k = j:length(X.element_dict), l = k:length(X.element_dict)
            if X.triplet_SRO[j, k, l, i] != 0
                if j == k && k == l && j == l
                    continue
                else
                    triplet_SRO_sum[i] += 1 - X.triplet_SRO[j, k, l, i]
                end
            end
        end
    end
    total_triplet_sum = 0
    for i = 1:length(triplet_SRO_sum)
        total_triplet_sum += 1/i*triplet_SRO_sum[i]
    end
    X.triplet_SRO = T
    return total_pair_sum, total_triplet_sum
end
