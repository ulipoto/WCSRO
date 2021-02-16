include("initialize.jl")
include("shuffle_config.jl")


"""
Calculates the pair- and triplet-SRO sum separately.
Inputs are a 'Backpack' and the sequence to be calculated.
"""
function sro(X::Backpack, sequence)
    A = copy(X.alpha)
    B = copy(X.beta)
    alpha_sum = 0
    beta_sum = 0
# Pair part: Going through all combinations (without double counts)
    for i = 1:length(sequence), j = i+1:length(sequence)
        if sequence[i][2] != 0 && sequence[j][2] != 0
            A[sequence[i][2], sequence[j][2], X.shell_matrix[i, j]] -=
            X.pair_prefactors[sequence[i][2], sequence[j][2], X.shell_matrix[i, j]]
        end
    end
# Summation (and substracting all pairs between to equal elements)
    for i = 1:size(A)[3]
        alpha_sum += 1/i*(sum(broadcast(abs, A[:, :, i])) - sum([abs(A[n, n, i]) for n in 1:size(A)[1]]))
    end

# Triplet part: This time, the tuple values of each element not only indicate
# the index, but also act as a key for the right mole fraction:
    for i = 1:length(sequence), j = 1:length(sequence), k = 1:length(sequence)
        if sequence[i][2] != 0 && sequence[j][2] != 0 && sequence[k][2] != 0
            si, sj, sk = X.triplet_shell_matrix[i, j, k]
            (si, sj, sk) == (0,0,0) && continue
            spa = sequence[i][2]
            spb = sequence[j][2]
            spc = sequence[k][2]
            xa = X.mole_fractions[sequence[i]]
            xb = X.mole_fractions[sequence[j]]
            xc = X.mole_fractions[sequence[k]]
            prefactor = X.trip_prefactors[(si, sj, sk)]/(xa*xb*xc)
            B[spa, spb, spc, si, sj, sk] -= prefactor
        end
    end
# Summation (and substracting all triplets between to equal elements)
    for i = 1:size(B)[4], j = 1:size(B)[5], k = 1:size(B)[6]
        beta_sum += 1/(3*i*j*k)*(sum(broadcast(abs, B[:, :, :, i, j, k])) - sum([abs(B[n, n, n, i, j, k]) for n in 1:size(B)[1]]))
    end
    return alpha_sum, beta_sum
end

# Helping structure for the results
struct SRO_results
    pairs::Vector
    triplets::Vector
    sro_sum::Vector
    pair_extrema::Vector{Tuple{Float64,Int64}}
    triplet_extrema::Vector{Tuple{Float64,Int64}}
    overall_extrema::Vector{Tuple{Float64,Int64}}
end

"""
Helps to evaluate the desired calculations and outputs the minima and maxima.
"""
function evaluate(X::Backpack, permutations, weight_factor::Float64)
    perms = permutations
    my_pairs = []
    my_triplets = []
    my_sros = []
    for i = 1:length(perms)
        current_sro = sro(X, perms[i])
        push!(my_pairs, current_sro[1])
        push!(my_triplets, current_sro[2])
        push!(my_sros, current_sro[1]+weight_factor*current_sro[2])
    end
    pair_extrema = [findmin(my_pairs), findmax(my_pairs)]
    triplet_extrema = [findmin(my_triplets), findmax(my_triplets)]
    overall_extrema = [findmin(my_sros), findmax(my_sros)]
    return SRO_results(my_pairs, my_triplets, my_sros, pair_extrema, triplet_extrema, overall_extrema)
end
