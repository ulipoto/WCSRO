using DelimitedFiles
using LinearAlgebra

include("structs.jl")
include("functions.jl")

function init_matrix(filename)
    #Necessary arrays for calculation
    sequence = []
    "The .xyz-file is being converted here:"
    total_atoms = parse(Int, readline(open(filename, "r")))
    dist_matrix = zeros(total_atoms, total_atoms)
    angle_matrix = zeros(total_atoms, total_atoms, total_atoms)
    u = (readdlm(filename,'\t', skipstart = 2, skipblanks = true))
    xyz = zeros(total_atoms, 3)
    dist_matrix = zeros(total_atoms, total_atoms)
    for i in 1:total_atoms
        l = split(u[i])
        push!(sequence, l[1])
        xyz[i, 1] = parse(Float64, l[2])
        xyz[i, 2] = parse(Float64, l[3])
        xyz[i, 3] = parse(Float64, l[4])
    end

    "All different elements are determined"
    elements, el_counter = element_counter(sequence)

    "Matrices for counting all possible pairs and triplet configurations
    are being created"
    #Pair matrix
    help_mat = []
    for i=1:length(elements)
        for j=i:length(elements)
            push!(help_mat, elements[i]*elements[j])
        end
    end
    pair_SRO = Dict(val => 0.0 for (idx, val) in enumerate(help_mat))
    #println(pair_SRO)

    #Triplet matrix
    help_mat = []
    for i=1:length(elements)
        for j=i:length(elements)
            for k=j:length(elements)
                push!(help_mat, elements[i]*elements[j]*elements[k])
            end
        end
    end
    triplet_SRO = Dict(val => 0.0 for (idx, val) in enumerate(help_mat))
    #println(triplet_SRO)

    """
    A distance matrix will be created. The value in the matrix at position [i,j]
    represents the norm of the distance between atom i and atom j. Furthermore,
    an array of all angles between the atoms i, j and k is being created, where
    index i is the atom that sits in the corner. Periodic boundary conditions
    (PBCs) have been taken into consideration. More on this in the pbc.jl docs.
    """
    dist_matrix, angle_matrix = pbc(xyz)
    X = Initializer(total_atoms, xyz, sequence, elements, el_counter, dist_matrix, angle_matrix, pair_SRO, triplet_SRO)
    return X
end
