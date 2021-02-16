include("structure.jl")
using LinearAlgebra

function initialize(file)
# PART OF YAML DATA EXTRACTION
    data = YAML.load_file(file)
    raw = data["structure"]
    s = Structure(
        hcat(raw["lattice"]...)',
        raw["species"],
        hcat(raw["fcoords"]...)',
        data["replicate"],
        data["composition"]
    )
    sc = supercell(s)
# Now, the shortest distance is being evualated by periodic boundary conditions
# and stored in 'vec_mat'
    rg = [1, 0, -1]
    shifts = (collect(Iterators.product(rg, rg, rg)))
    shifts = reshape(shifts, :, 1)
    cart_shifts = []
    distance_matrix = zeros(size(sc.fcoords)[1], size(sc.fcoords)[1])
    vec_mat = zeros(size(sc.fcoords)[1], size(sc.fcoords)[1], 3)
    min_dist = 1e-6
    for i in shifts
        push!(cart_shifts, collect(i))
    end

    for i in 1:size(sc.fcoords)[1]
        for j in 1:size(sc.fcoords)[1]
            if i == j
                continue
            else
            d_vecs = []
            for sh = 1:length(cart_shifts)
                push!(d_vecs, sc.lattice*reshape(sc.fcoords[j, :] + cart_shifts[sh] - sc.fcoords[i, :], (3, 1)))
            end
            for k = 1:length(d_vecs)
                for l = 1:length(d_vecs[k])
                    d_vecs[k][l] = round(d_vecs[k][l], digits = 5)
                end
            end
            min_vec = [Inf, Inf, Inf]
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
        for j = 1:size(vec_mat)[1]
            distance_matrix[i, j] = round(norm(vec_mat[i, j, :]), digits = 4) #changed round(digits 3!)
        end
    end
# Structures from the functions below are created here.
    distances = (sort!(unique(distance_matrix)))
    distance_dict = Dict(dist => index-1 for (index, dist) in enumerate(distances))
    shell_matrix = map(num -> distance_dict[num], distance_matrix)
    shells = coord_shells(shell_matrix)
    triplet_shell_matrix = triplet_shell(shell_matrix)
    map!(t -> cutoff(t), triplet_shell_matrix, triplet_shell_matrix) #induces cutoff
    mole_fracs, elem_dict = mole_fractions(sc)
    pair_prefactors = calc_pair_pref(mole_fracs, shells, shell_matrix)
    trip_prefactors = calc_trip_prefs(triplet_shell_matrix, shells)
# Finally, the two arrays for SRO value storage
    alpha = ones(length(mole_fracs), length(mole_fracs), length(shells))
    beta = ones(length(mole_fracs), length(mole_fracs), length(mole_fracs), length(shells), length(shells), length(shells))
# DONE :-)
    return Backpack(sc, elem_dict, mole_fracs, shell_matrix, triplet_shell_matrix, alpha, beta, pair_prefactors, trip_prefactors)
end

function cutoff((s1, s2, s3), cutoff=30)
# Simple cutoff help function
    any([s1 >= cutoff, s2 >= cutoff, s3 >= cutoff ]) && return (0,0,0)
    any([s1 == 0, s2 == 0, s3 == 0]) && return (0,0,0)
    return (s1, s2, s3)
end

function triplet_shell(shell_matrix::Array, cutoff = Inf)
# constructs a 3D-matrix containing a tuple of all three shells between
# atoms i, j and k
    sm = shell_matrix
    triplet_shell_matrix = [(sm[i, j], sm[j, k], sm[i, k]) for i=1:size(sm)[1], j=1:size(sm)[1], k=1:size(sm)[1]]
    return triplet_shell_matrix
end

function coord_shells(shell_matrix::Array)
# returns a dict with the shell number as key and the number of appearing atoms
# (with pbc) as values
    shell_dict = Dict(i => 0.0 for i = 1:maximum(shell_matrix))
    for i = 1:size(shell_matrix)[1], j = i:size(shell_matrix)[1]
        shell_matrix[i, j] == 0 && continue
        shell_dict[shell_matrix[i, j]] += 2/size(shell_matrix)[1]
    end
    return shell_dict
end

function calc_trip_prefs(triplet_shell_matrix::Array, M::Dict)
# returns a dictionary with all occurring triplets as keys and their corresponding
# prefactors
    tt = unique(triplet_shell_matrix) #all triplet types
    triplet_prefactors = Dict(map(t-> trip_pref_helper(t, M, size(triplet_shell_matrix)[1]), tt))
    return triplet_prefactors
end

function trip_pref_helper(s::Tuple, M::Dict, sc_size)
# Calculates the fixed part ν and rules out cases with a zero
    if any(i -> i == 0, s) == true
        return (0,0,0) => nothing
    else
        return s => 1/(3*sc_size)*(1/(M[s[1]]*M[s[2]])+1/(M[s[1]]*M[s[3]])+1/(M[s[2]]*M[s[3]]))
    end
end

function calc_pair_pref(mole_fractions, shells, shell_matrix)
# returns an 3D-array with all prefactors. Prefactor for a pair between elements
# 1 and 3 - 4 shells apart- is found at pair_prefactors[1, 3, 4] (and [3, 1, 4])
    pair_prefactors = zeros(length(mole_fractions), length(mole_fractions), length(shells))
    for i = 1:length(shells)
        for (idx1, val1) in enumerate(mole_fractions)
            if val1[1][2] != 0
                for (idx2, val2) in enumerate(mole_fractions)
                    if val2[1][2] != 0
                        #Factor 2 because there is no double counting later on
                        pair_prefactors[val1[1][2], val2[1][2], i] = 2/(val1[2]*val2[2]*shells[i]*size(shell_matrix)[1])
                    end
                end
            end
        end
    end
    return pair_prefactors
end

function mole_fractions(s::Structure)
# Returns mole fraction and an dictionary with all elements in it. The values
# in element_dict are the keys for mole_fractions
    mf = Dict()
    elem_dict = Dict()
    seq = s.species
    elem_num = 1
    for (idx1, val1) in enumerate(s.composition)
        sub_sum = (sum(values(val1[2])))
        for (idx2, val2) in enumerate(val1[2])
            if val2[1] == "Vac" #Vacancy exception: returns 0
                mf[(idx1, 0)] = val2[2]/sub_sum
                elem_dict[val2[1]] = (idx1, 0)
                i = 1
                j = 1
            else
                mf[(idx1, elem_num)] = val2[2]/sub_sum #key is sublattice + element number
                elem_dict[val2[1]] = (idx1, elem_num) #key above is the value here
                i = 1
                j = 1
                elem_num += 1
            end
            for i = 1:length(seq)
                if seq[i] == val1[1] && j <= val2[2]
                    seq[i] = val2[1]
                    j += 1
                end
            end
            i += 1
        end
    end
    return mf, elem_dict
end
