using IterTools
using LinearAlgebra

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
    for i = 1:length(elements)
        element_dict[elements[i]] = element_counter[i]
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

    """
    The shortest distances between all atom pairs are calculated and stored in
    'distance matrix' by usage of periodic boundariy conditions (PBC).
    """
    rg = [1, 0, -1]
    shifts = (collect(Iterators.product(rg, rg, rg)))
    shifts = reshape(shifts, :, 1)
    cart_shifts = []
    distance_matrix = zeros(length(fcoords), length(fcoords))
    vec_mat = zeros(length(fcoords), length(fcoords), 3)
    min_dist = 1e-6
    for i in shifts
        push!(cart_shifts, lattice * collect(i))
    end
    for (i, fc1) in enumerate(fcoords)
        for (j, fc2) in enumerate(fcoords)
            if i == j
                continue
            else
            cc1 = lattice * fc1
            cc2 = lattice * fc2
            d_vecs = []
            d_vecs = [cc2 + sh - cc1 for sh in cart_shifts]
            min_vec = [1000, 1000, 1000]
                for k in d_vecs
                    if norm(k) < norm(min_vec) && norm(k) > min_dist
                        vec_mat[i, j, :] = k
                    end
                end
            end
        end
    end
    for i = 1:size(vec_mat)[1]
        for j = 1:size(vec_mat)[1]
            distance_matrix[i, j] = norm(vec_mat[i, j, :])
        end
    end

    "Matrices for counting all possible pairs and triplet configurations
    are being created"
    #Pair matrix
    pair_SRO = Dict()
    for i=1:length(elements)
        for j=i:length(elements)
            pair_SRO[elements[i], elements[j]] = 0
        end
    end

    #Triplet matrix
    triplet_SRO = Dict()
    for i=1:length(elements)
        for j=i:length(elements)
            for k=j:length(elements)
                triplet_SRO[elements[i], elements[j], elements[k]] = 0
            end
        end
    end
    X = Initializer(fcoords, lattice, element_dict, distance_matrix, pair_SRO, triplet_SRO)
    return X
end

"""
get_pair_key(D::Dict, a)

The function returns the corresponding key in the SRO dictionary for
a given element pair.
"""
function get_pair_key(D::Dict, a::Tuple)
    perms = collect(permutations(a, 2))
    for i = 1:length(perms)
        if haskey(D, (perms[i][1], perms[i][2])) == true
            return (perms[i][1], perms[i][2])
        end
    end
end

"""
get_pair_key(D::Dict, a)

The function returns the corresponding key in the SRO dictionary for
a given element triplet.
"""
function get_triplet_key(D::Dict, a::Tuple)
    perms = collect(permutations(a, 3))
    for i = 1:length(perms)
        if haskey(D, (perms[i][1], perms[i][2], perms[i][3])) == true
            return (perms[i][1], perms[i][2], perms[i][3])
        end
    end
end
