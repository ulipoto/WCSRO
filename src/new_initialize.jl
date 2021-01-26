using IterTools
using LinearAlgebra
include("structs.jl")
include("sro_calculator.jl")

file = "POSCAR.mp-88_AB" # Please note, that the file has to be in the same folder.

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
    filename = file
    element_counter = [] #list of occurrence of all elements
    elements = [] #list of all different occurring elements
    element_dict = Dict() # Dictionary with the element as a key and a tuple of
                          # an assigned index and the occurrence of the element.
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

    # Fractional coordinates and lattice parameters are extracted
    for (num, line) in enumerate(eachline(filename))
        if num > 8
            push!(fcoords, [parse(Float64, match(fcoords_regex, line)[1]),  parse(Float64, match(fcoords_regex, line)[2]),
            parse(Float64, match(fcoords_regex, line)[3])])
        end
    end
    for (num, line) in enumerate(eachline(filename))
        if num > 2 && num < 6
            push!(lattice, parse(Float64, match(fcoords_regex, line)[1]), parse(Float64, match(fcoords_regex, line)[2]),
            parse(Float64, match(fcoords_regex, line)[3]))
        end
    end
    lattice = reshape(lattice, 3, 3)
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
            d_vecs = zeros(length(cart_shifts), 2)
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
    fcc_shells = [12, 6, 24, 12, 24, 8, 48] #Coordination numbers for fcc: Max = 134!

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

    "Weight matrices are being create. This should help later on, when a lot of
    iterations are necessary. The pair weights are stored in an N x N x S matrix,
    where N is the number of different elements in the system and S the number of
    coordinations shells. Thus, in a binary system, a pair between element A and B,
    where B is in the third corrdination shell, is stored at position [1, 2, 3].
    Same is true for triplets. There, however, the number of the two different
    coordination shells is added together."
    for i = 1:length(elements)
        element_dict[elements[i]] = (i, element_counter[i])
    end

    pair_weight_matrix = zeros(length(element_dict), length(element_dict), length(distance_dict))
    for i = 1:length(distance_dict)
        for (idx1, val1) in enumerate(element_dict)
            for (idx2, val2) in enumerate(element_dict)
                pair_weight_matrix[val1[2][1], val2[2][1], i] = 1/(length(fcoords)*fcc_shells[i]*(val1[2][2]/length(fcoords))*(val2[2][2]/length(fcoords)))
            end
        end
    end
    max_shells = 4
    triplet_weight_matrix = zeros(length(element_dict), length(element_dict), length(element_dict), max_shells)
    for i = 1:3, j = 1:3, k = 1:3
        for (idx1, val1) in enumerate(element_dict)
            for (idx2, val2) in enumerate(element_dict)
                for (idx3, val3) in enumerate(element_dict)
                    if i + j + k - 2 <= max_shells
                        if i == j
                            triplet_weight_matrix[val1[2][1], val2[2][1], val3[2][1], i + j + k - 2] = 1/(length(fcoords)*fcc_shells[i]*(fcc_shells[j]-1)*
                            (val1[2][2]/length(fcoords))*(val2[2][2]/length(fcoords))*(val3[2][2]/length(fcoords)))
                        else
                            triplet_weight_matrix[val1[2][1], val2[2][1], val3[2][1], i + j + k - 2] = 1/(length(fcoords)*fcc_shells[i]*fcc_shells[j]*
                            (val1[2][2]/length(fcoords))*(val2[2][2]/length(fcoords))*(val3[2][2]/length(fcoords)))
                        end
                    end
                end
            end
        end
    end

    "Matrices for counting all possible pairs and triplet configurations
    are being created"
    pair_SRO = zeros(length(element_dict), length(element_dict), length(distance_dict))
    triplet_SRO = zeros(length(element_dict), length(element_dict), length(element_dict), max_shells)
    X = Initializer(element_dict, distances, shell_matrix, pair_weight_matrix, triplet_weight_matrix, pair_SRO, triplet_SRO)
    return X
end
x = initialize(file)
