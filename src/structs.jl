"""
The Struct 'Backpack' is a struct, which shall contain all variables that
have to be calculated just once. It serves to get rid of potentially disturbing
global variables.
"""

mutable struct Backpack
    element_dict::Dict ## Dictionary with the element as a key and a tuple of
                       ## an assigned index and the occurrence of the element.
    distances::Array ## All occurring distances
    all_pairs::Array ## An array with all possible pairs
    all_triplets::Array ## An array with all possible triplets
    pair_weight_matrix::Array ## Weights for all pairs in different shells
    triplet_weight_matrix::Array ## Weights for all triplets in different shells
    pair_SRO::Array ## Dictionary for all pair interactions
    triplet_SRO::Array ## Dictionary for all pair interactions
end
