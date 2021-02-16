include("structure.jl")
include("initialize.jl")
include("shuffle_config.jl")
include("sro_calc.jl")
include("plotting.jl")

"""
WELCOME TO MY SRO CALCULATOR.
This script calculates Warren-Cowley short range order parameters for pairs
and triplets in arbitrary alloys. The code works for an arbitrary number of
different elements. This file shall guide you through the necessary steps.
"""

# INPUT
# Please declare your input YAML-file:

file = "testfile.yaml"

# This will create all structures for a SRO calculation and stored it in an
# object of type 'Backpack'.

my_backpack = initialize(file)

# Now, the SRO parameters are being calculated while the atomic configuration is
# shuffled in between. The input number below determines the number of
# non-repetitive shuffles:

shuffles = 100
perms = shuffle_atoms(my_backpack, shuffles)
#xyz_animation(my_backpack, unique_perms[2])

# Alternatively, if you want to generate aLL permutations, you can use this function
# below. You find the documentations for this functions in 'shuffle_sonfig.jl' as well.
# Please keep in mind that this can be very CPU-intensive and is therefore not
# suitable for too large systems. But feel free to try it out. ;-)
# @time begin
    # perms = all_permutations(sequence(my_backpack))
# end


# Now, we start SRO calculation. With the function 'evaluate', SROs are calculated for
# all created permutations and the minima and maxima for pairs, triplets and a
# sum of both are returned separately. The value in the function is the triplet weight factor,
# which defines the contribution of the triplets
results = evaluate(my_backpack, perms, 0.2)

# Finally, see what you did! :-)
# You can either plot some reults, like here the pairs...
plot_configs(results.pairs)
#...or you plot the whole structure, for example the minimum for triplets...
#xyzplot(my_backpack, perms[results.triplet_extrema[1][2]], 15.0)
#...you can even create a little gif, if you want:
#xyz_animation(my_backpack, perms[results.triplet_extrema[1][2]])
