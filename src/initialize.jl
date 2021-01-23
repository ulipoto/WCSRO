using IterTools
using LinearAlgebra
using YAML
include("structs.jl")

"""
This file contains a function 'initmatrix' for processing an input file.
The input file can either be of type .yaml or POSCAR.
"""

"""
initialize(file)

Necessary structures for SRO calculation are created and stored in a object
of type "Backpack".
"""
function initialize(file)

# PART OF POSCAR DATA EXTRACTION
    filename = file
    if occursin("POSCAR", filename) == true
        element_counter = [] #list of occurrence of all elements
        elements = [] #list of all different occurring elements
        element_dict = Dict() # Dictionary of all elements and their respective occurrence.
        fcoords = [] #All fractional coordinates
        lattice = [] #lattice parameter array
        s = read(open(filename, "r"), String)
        # Regex to find the number of occurrence for a given element
        atomic_int_regex = r"(?<=\s|^)\d+(?=\s|$)"
        element_regex = r"\b[^\d\W]+\b"
        fcoords_regex = r"(?P<fx>\d+\.\d+)\s+(?P<fy>\d+\.\d+)\s+(?P<fz>\d+\.\d+).*"

        # Elements and their respective occurrence are extracted from the file.
        # A dictionaray will be created subsequently.
        for i in eachmatch(atomic_int_regex, s)
            push!(element_counter, parse(Int64, i.match))
        end
        for i in eachmatch(element_regex, s)
            if i.match == "direct"
                break
            else
                push!(elements, i.match)
            end
        end
        for (idx, val) in enumerate(elements)
            element_dict[val[1]] = (idx, element_counter[idx])
        end
        # Fractional coordinates and lattice parameters are extracted
        for (num, line) in enumerate(eachline(filename))
            if num > 8
                push!(fcoords, [parse(Float64, match(fcoords_regex, line)[1]),  parse(Float64, match(fcoords_regex, line)[2]),
                parse(Float64, match(fcoords_regex, line)[3])])
            end
        end
        for (num, line) in enumerate(eachline(filename))
            if num > 2 && num < 6
                push!(lattice, parse(Float64, match(fcoords_regex, line)[1]),  parse(Float64, match(fcoords_regex, line)[2]),
                parse(Float64, match(fcoords_regex, line)[3]))
            end
        end
        lattice = reshape(lattice, 3, 3)
    else
    # PART OF YAML DATA EXTRACTION
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
    end

# Now, the shortest distance is being evualated by periodic boundary conditions
# and stored in 'vec_mat'
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

# The distance matrix is being created out of the shortest distances
    for i = 1:size(vec_mat)[1]
        for j = i:size(vec_mat)[1]
            distance_matrix[i, j] = norm(vec_mat[i, j, :])
        end
    end
# Now, all different distances are evaluated and a corresponding shell matrix
# and a matrix containing all triplets is made.
# First for pairs....
    distances = zeros(1)
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
    pair_shell_matrix = zeros(Int64, size(distance_matrix)[1], size(distance_matrix)[2])
    for i = 1:size(distance_matrix)[1], j = i:size(distance_matrix)[2]
        pair_shell_matrix[i, j] = Int(distance_dict[distance_matrix[i, j]])
    end
    all_pairs = []
    exp_pairs = zeros(length(distance_dict))
    for i = 1:size(pair_shell_matrix)[1]
        for j = 1:size(pair_shell_matrix)[2]
            if pair_shell_matrix[i, j] != 0
                exp_pairs[pair_shell_matrix[i, j]] += 1
                push!(all_pairs, [i, j, pair_shell_matrix[i, j]])
            end
        end
    end
    while exp_pairs[end] == 0
        pop!(exp_pairs)
    end
# ...now for triplets
    all_triplets = []
    exp_triplets = zeros(3)
    for i = 1:size(pair_shell_matrix)[1]
        for j = 1:size(pair_shell_matrix)[2]
            if pair_shell_matrix[i, j] > 0 && pair_shell_matrix[i, j] < 3
                for k = j:size(pair_shell_matrix)[1]
                    if pair_shell_matrix[i, k] > 0 && pair_shell_matrix[i, k] < 3
                        exp_triplets[pair_shell_matrix[i, j] + pair_shell_matrix[i, k] - 1] += 1
                        push!(all_triplets, [i, j, k, pair_shell_matrix[i, j] + pair_shell_matrix[i, k] - 1])
                    end
                end
            end
        end
    end

"Weight matrices are being created. This should help at SRO calculation.
The pair weights are stored in an N x N x S matrix, where N is the number of
different elements in the system and S the number of coordinations shells.
Thus, in a binary system, a pair between element A (=1) and B (=2), where B is
in the third corrdination shell, is stored at position [1, 2, 3].
Same is true for triplets. There, however, there are only three shells: First
with both in the first shell, second with first/second and third with
second/second."
    pair_weight_matrix = zeros(length(element_dict), length(element_dict), length(exp_pairs))
    for i = 1:length(exp_pairs)
        for (idx1, val1) in enumerate(element_dict)
            for (idx2, val2) in enumerate(element_dict)
                if val1[2][1] == val2[2][1]
                    pair_weight_matrix[val1[2][1], val2[2][1], i] = 1/(exp_pairs[i]*(val1[2][2]/length(fcoords))*(val1[2][2]/length(fcoords)))
                else
                    pair_weight_matrix[val1[2][1], val2[2][1], i] = 1/(2*exp_pairs[i]*(val1[2][2]/length(fcoords))*(val1[2][2]/length(fcoords)))
                end
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
                        triplet_weight_matrix[vals[1], vals[2], vals[3], i] = 1/(exp_triplets[i]*(val1[2][2]/length(fcoords))*
                        (val2[2][2]/length(fcoords))*(val3[2][2]/length(fcoords)))
                    elseif vals[1] != vals[2] && vals[1] != vals[3] && vals[2] != vals[3]
                        triplet_weight_matrix[vals[1], vals[2], vals[3], i] = 1/(6*exp_triplets[i]*(val1[2][2]/length(fcoords))*
                        (val2[2][2]/length(fcoords))*(val3[2][2]/length(fcoords)))
                    else
                        triplet_weight_matrix[vals[1], vals[2], vals[3], i] = 1/(3*exp_triplets[i]*(val1[2][2]/length(fcoords))*
                        (val2[2][2]/length(fcoords))*(val3[2][2]/length(fcoords)))
                    end
                end
            end
        end
    end
"Matrices for counting all possible pairs and triplet configurations
are being created"
    pair_SRO = zeros(length(element_dict), length(element_dict), length(distance_dict))
    triplet_SRO = zeros(length(element_dict), length(element_dict), length(element_dict), length(exp_triplets))
    X = Backpack(element_dict, distances, all_pairs, all_triplets, pair_weight_matrix, triplet_weight_matrix, pair_SRO, triplet_SRO)
    return X
end

#file = "ab.yml"
file = "POSCAR.mp-88_AB"
X = initialize(file)
X
