module WCSRO

include("structs.jl")
include("initialize.jl")
include("atom_shuffle.jl")
include("sro_calculator.jl")
include("plotting.jl")

"""
WELCOME TO MY SRO CALCULATOR.
This script calculates Warren-Cowley short range order parameters for pairs
and triplets in arbitrary alloys. The code works for an arbitrary number of
different elements. This file shall guide you through the necessary steps.
"""

# INPUT
file = "POSCAR.mp-1184260_Fe3Cu" # Please note, that the file has to be in the same folder.

# This will create all necessary arrays and dictionaries for SRO calculations.
# All values are stored in a struct. Here, this struct is called 'X'.
X = initialize(file)
# Now, the SRO parameters are being calculated while the atomic configuration
# is shuffled in between. The input number below determines the number
# of non-repetitive shuffles.
shuffles = 50
unique_perms = shuffle_atoms(X, shuffles)
best_pairs, best_triplets = optimize(X, unique_perms)
println(best_pairs)
println(best_triplets)
# Finally, the best pair structure from our sample size is being plotted
xyzplot(X, unique_perms[1, :])
end
