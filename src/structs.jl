"""
The Struct 'Initializer' is a struct, which shall contain all variables that
have to be calculated just once. It serves to get rid of potentially disturbing
global variables.
"""

struct Initializer
    fcoords::Array ## Array with all fractional coordinates coordinates
    lattice::Array ## Array with lattice constants
    element_dict::Dict ## Dictionary for the number of atoms
    distance_matrix::Array ## Array with shortest distance between i and j (with PBC)
    pair_SRO::Dict ## Dictionary for all pair interactions
    triplet_SRO::Dict ## Dictionary for all pair interactions
end
