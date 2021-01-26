using IterTools
using LinearAlgebra
using YAML
include("structs.jl")

"""
This file contains a function 'initmatrix' for processing an input file.
Furthermore, there are two helping functions located below for creating correct
dictionary keys. These functions are necessary in the file 'sro_calculator.jl'.
"""

"""
initialize(file)

Necessary arrays and regex definitions are created in this function. The imported
file will be opened, informations like elements and coordinates are extracted and
distributed to corresponding arrays.
"""
function initialize(file)
    data = YAML.load_file(file)
    lattice = data["lattice parameters"]
    comp = data["composition"]
    rep = data["replication"]
    element_dict = Dict()
    fcoords = []

    while length(lattice) < 3
        push!(lattice, lattice[1])
    end
    lattice = Diagonal(lattice)
    uc = [0 0 0; 0.5 0.5 0; 0.5 0 0.5; 0 0.5 0.5]
    uc = uc/maximum(rep)
    xtrans = collect(0:1/rep[1]:0.9999)
    ytrans = collect(0:1/rep[2]:0.9999)
    ztrans = collect(0:1/rep[3]:0.9999)
    for i = 1:length(xtrans), j = 1:length(ytrans), k = 1:length(ztrans)
        for m = 1:size(uc)[1]
            push!(fcoords, [uc[m, 1]+xtrans[i], uc[m, 2]+xtrans[j], uc[m, 3]+xtrans[k]])
        end
    end
    comp_vec = []
    for (idx, val) in enumerate(comp)
        x = round(val[2]/100*length(fcoords), digits = 0)
        push!(comp_vec, Int(x))
    end
    if sum(comp_vec) > length(fcoords)
        comp_vec[end] -= 1
    elseif sum(comp_vec) < length(fcoords)
        comp_vec[end] += 1
    end
    for (idx, val) in enumerate(comp)
        element_dict[val[1]] = (idx, comp_vec[idx])
    end

    rg = [1, 0, -1]
    shifts = (collect(Iterators.product(rg, rg, rg)))
    shifts = reshape(shifts, :, 1)
    cart_shifts = []
    distance_matrix = zeros(length(fcoords), length(fcoords))
    vec_mat = zeros(length(fcoords), length(fcoords), 3)
    min_dist = 1e-6
    for i in shifts
        push!(cart_shifts, collect(i))
    end
    for (i, fc1) in enumerate(fcoords)
        for (j, fc2) in enumerate(fcoords)
            if i == j
                continue
            else
            d_vecs = []
            d_vecs = [fc2 + sh - fc1 for sh in cart_shifts]
            for k = 1:length(d_vecs)
                for l = 1:length(d_vecs[k])
                    d_vecs[k][l] = round(d_vecs[k][l], digits = 5)
                end
            end
            min_vec = [1000, 1000, 1000]
                for k in d_vecs
                    if norm(k) < norm(min_vec) && norm(k) > min_dist
                        vec_mat[i, j, :] = k
                        min_vec = k
                    end
                end
            end
        end
    end
    all_entries = zeros(1)
    for i = 1:size(vec_mat)[1]
        for j = 1:size(vec_mat)[2]
            for m = 1:length(all_entries)
                if isapprox(norm(vec_mat[i, j, :]), all_entries[m], atol = 1e-3) == true
                    break
                end
                if m == length(all_entries)
                    push!(all_entries, norm(vec_mat[i, j, :]))
                end
            end
        end
    end
    sort!(all_entries)
    #The distance matrix is being created
    for i = 1:size(vec_mat)[1]
        for j = 1:size(vec_mat)[1]
            distance_matrix[i, j] = norm(vec_mat[i, j, :])
        end
    end
    # Now, all different distances are evaluated and a corresponding shell matrix
    # is made.
    distances = zeros(1)
    distance_dict = Dict()
    fcc_shells = zeros(6) #Coordination numbers for fcc: Max = 134!

    for i = 1:size(distance_matrix)[1], j = 1:size(distance_matrix)[1]
        if any(x -> isapprox(x, distance_matrix[i, j], atol = 1e-2), distances) == false
            push!(distances, distance_matrix[i, j])
        end
    end
    sort!(distances)
    distance_dict = Dict()
    for i = 1:length(distances)
        distance_dict[distances[i]] = i-1
    end
    shell_matrix = zeros(Int64, size(distance_matrix)[1], size(distance_matrix)[2])
    for i = 1:size(distance_matrix)[1], j = i:size(distance_matrix)[2]
        shell_matrix[i, j] = Int(distance_dict[distance_matrix[i, j]])
    end

    fcc_shells = [12, 6, 24, 12, 24, 8, 48]
    exp_triplets = zeros(3)
    for i = 1:size(shell_matrix)[1]
        for j = 1:size(shell_matrix)[2]
            if shell_matrix[i, j] > 0 && shell_matrix[i, j] < 3
                for k = j:size(shell_matrix)[1]
                    if shell_matrix[i, k] > 0 && shell_matrix[i, k] < 3
                        exp_triplets[shell_matrix[i, j] + shell_matrix[i, k] - 1] += 1
                    end
                end
            end
        end
    end
    "Weight matrices are being created. This should help later on, when a lot of
    iterations are necessary. The pair weights are stored in an N x N x S matrix,
    where N is the number of different elements in the system and S the number of
    coordinations shells. Thus, in a binary system, a pair between element A (=1)
    and B (=2), where B is in the third corrdination shell, is stored at position
    [1, 2, 3].
    Same is true for triplets. There, however, the number of the two different
    coordination shells is added together."
    pair_weight_matrix = zeros(length(element_dict), length(element_dict), length(distance_dict))
    for i = 1:length(distance_dict)
        for (idx1, val1) in enumerate(element_dict)
            for (idx2, val2) in enumerate(element_dict)
                pair_weight_matrix[val1[2][1], val2[2][1], i] = 1/(length(fcoords)*fcc_shells[i]*(val1[2][2]/length(fcoords))*(val2[2][2]/length(fcoords)))
            end
        end
    end
    
    triplet_weight_matrix = zeros(length(element_dict), length(element_dict), length(element_dict), length(exp_triplets))
    for i = 1:length(exp_triplets)
        for (idx1, val1) in enumerate(element_dict)
            for (idx2, val2) in enumerate(element_dict)
                for (idx3, val3) in enumerate(element_dict)
                    vals = [val1[2][1], val2[2][1], val3[2][1]]
                    sort!(vals)
                    if vals[1] == vals[2] && vals[1] == vals[3] && vals[2] == vals[3]
                        triplet_weight_matrix[vals[1], vals[2], vals[3], i] = 1/(exp_triplets[i]/(length(element_dict)^3))
                    elseif vals[1] != vals[2] && vals[1] != vals[3] && vals[2] != vals[3]
                        triplet_weight_matrix[vals[1], vals[2], vals[3], i] = 1/(6*exp_triplets[i]/(length(element_dict)^3))
                    else
                        triplet_weight_matrix[vals[1], vals[2], vals[3], i] = 1/(3*exp_triplets[i]/(length(element_dict)^3))
                    end
                end
            end
        end
    end
    "Matrices for counting all possible pairs and triplet configurations
    are being created"
    pair_SRO = zeros(length(element_dict), length(element_dict), length(distance_dict))
    triplet_SRO = zeros(length(element_dict), length(element_dict), length(element_dict), length(exp_triplets))
    X = Initializer(element_dict, distances, shell_matrix, pair_weight_matrix, triplet_weight_matrix, pair_SRO, triplet_SRO)
    return X
end

file = "ab.yml" # Please note, that the file has to be in the same folder.

X = initialize(file)
println("HI")
