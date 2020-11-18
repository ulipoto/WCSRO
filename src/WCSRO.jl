module WCSRO

include("init.jl")
include("functions.jl")
include("structs.jl")
include("atom_shuffle.jl")
include("sro_calculator.jl")
include("optimizer.jl")
include("plotting.jl")
"""
WELCOME TO MY SRO CALCULATOR.
This script calculates Warren-Cowley short range order parameters for pairs
and triplets in alloys. The code works for an arbitrary number of different
elements. This file shall guide you through the necessary steps.
"""

# INPUT
file = "Al3TiNi.xyz" # Please note, that the file has to be in the same folder.

# This will create all necessary arrays and dictionaries for SRO calculations.
# All values are stored in a struct. Here, I call it 'X'.
X = init_matrix(file)

# Now, we calculate the SRO parameters while shuffling our atomic configuration
# in between. The input number (integer!) below determines the number of non-
# repetitive shuffles.
shuffles = 50
my_permutations = shuffle_atoms(X.sequence, shuffles)
best_pairs, best_triplets = optimize(X, my_permutations)

# Finally, the best pair structure from our sample size is being plotted
xyzplot(X, best_pairs)
end
