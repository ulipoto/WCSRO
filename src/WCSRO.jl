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
file = "POSCAR.mp-865235_Ti6Al16Ni7" # Please note, that the file has to be in the same folder.

"""
This will create all necessary arrays and dictionaries for SRO calculations.
All values are stored in a struct. Here, this struct is called 'X'.
"""
X = initialize(file)

"""
Now, the SRO parameters are being calculated while the atomic configuration is
shuffled in between. The input number below determines the number of
non-repetitive shuffles:
"""
shuffles = 100
unique_perms = shuffle_atoms(X, shuffles)

"""
Alternatively, if you want to generate ALL permutations, you can use this function
below. You find the documentations for this functions in 'atom_shuffle.jl' as well.
Please keep in mind that this can be very CPU-intensive and is therefore not suit-
-able for too large systems. But feel free to try it out. ;-)
"""
#all_perms = all_permutations(sequence(X))

"""
Now, we create all corresponding SRO values, which are simple sums of pairs and
triplets without any weighing factors. So please keep in mind that these values
are only qualitative measurements for order.
"""
sro_list = []
for i = 1:size(unique_perms)[1]
    push!(sro_list, sro(X, unique_perms[i, :]))
end

# Finally, the best pair structure from our sample size is being plotted
xyzplot(X, unique_perms[findmin(sro_list)[2], :])
end
