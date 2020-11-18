"""The Struct 'Initializer' is a struct, which shall contain all variables, that
have to be calculated just once."
"""

struct Initializer
    total_atoms::Int ## The number of total atoms
    xyz::Array ## Array with cartesian coordinates
    sequence::Array ## The sequence of atoms from the input file
    elements::Array ## A list of all occurring elements
    element_dict::Dict ## Dictionary for the number of atoms
    distance_matrix::Array ## Array with shortest distance between i and j (with PBC)
    angles::Array ## Array with angles between i, j and k (with PBC)
    pair_SRO::Dict ## Dictionary for all pair interactions
    triplet_SRO::Dict ## Dictionary for all pair interactions
end
