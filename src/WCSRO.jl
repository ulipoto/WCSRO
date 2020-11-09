module WCSRO

include("init.jl")
include("atom_shuffle.jl")
include("sro_calculator.jl")
include("plotting.jl")

"Welcome to my WCSRO calculator!
At the moment, most of the parameters are already chosen. In the future, the
parameters should either be an input via an UI or changeable here in this file."

#INPUT
file = "Al3Ti.xyz" #Please note, that the file has to be in the same folder.

#This will create all necessary arrays for SRO calculation.
init_matrix(file)

#The input determines the number of non-repetitive shuffles of the different atoms
iterations = 50

shuffle_atoms(iterations)
SRO(unique_perms, iterations)
println(SRO_array)

#Plotting the current best structure
xyzplot(unique_perms[findmin(SRO_array)[2]])

"For now, it just prints out the SRO parameters for a binary alloy. The ternary
part is already included and fundamentally tested, however the example file is
a binary alloy anyway. As I already said, this file needs further improvements,
it just should serve as a little proof of concept :-)"
end
