"""
The Struct 'Initializer' is a struct, which shall contain all variables that
have to be calculated just once. It serves to get rid of potentially disturbing
global variables.
"""

struct Initializer
    element_dict::Dict ## Dictionary with the element as a key and a tuple of
                       ## an assigned index and the occurrence of the element.
    distances::Array ## All occurring distances
    shell_matrix::Array ## Distance matrix expressed by the shell
    pair_weight_matrix::Array ## Weights for all pairs in different shells
    triplet_weight_matrix::Array ## Weights for all triplets in different shells
    pair_SRO::Array ## Dictionary for all pair interactions
    triplet_SRO::Array ## Dictionary for all pair interactions
end
