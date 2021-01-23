module WCSRO
include("structs.jl")
include("initialize.jl")
include("atom_shuffle.jl")
include("sro_calc.jl")
include("plotting.jl")

"""
WELCOME TO MY SRO CALCULATOR.
This script calculates Warren-Cowley short range order parameters for pairs
and triplets in arbitrary alloys. The code works for an arbitrary number of
different elements. This file shall guide you through the necessary steps.
"""

# INPUT
# The input can either be of type '.yaml' or a POSCAR-file:
file = "ab.yml"

# This will create all structures for a SRO calculation and stored it in an
# object of type 'Backpack'.
my_backpack = initialize(file)


#Now, the SRO parameters are being calculated while the atomic configuration is
#shuffled in between. The input number below determines the number of
#non-repetitive shuffles:
shuffles = 100
unique_perms = shuffle_atoms(my_backpack, shuffles)
@time shuffle_atoms(my_backpack, shuffles)
unique_perms


#Alternatively, if you want to generate ALL permutations, you can use this function
#below. You find the documentations for this functions in 'atom_shuffle.jl' as well.
#Please keep in mind that this can be very CPU-intensive and is therefore not
#suitable for too large systems. But feel free to try it out. ;-)
# @time begin
    # all_perms = all_permutations(sequence(my_backpack))
# end
# println(size(all_perms))


#Now, we create all corresponding SRO values, which are simple sums of pairs and
#triplets without any weighing factors. So please keep in mind that these values
#are only qualitative measurements for order.

sro_list = []
for i = 1:size(unique_perms)[1]
    push!(sro_list, @time sro(my_backpack, unique_perms[i, :]))
end
end
