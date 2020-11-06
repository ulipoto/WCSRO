module WCSRO
include("init.jl")
include("sro_calculator.jl")
include("atom_shuffle.jl")


"Welcome to my WCSRO calculator!
At the moment, most of the parameters are already chosen. In the future, the
parameters should either be an input via an UI or changeable here in this file."

#This will create all necessary arrays for SRO calclulation.
init_matrix()

#The input determines the number of non-repetitive shuffles of the different atoms
shuffle_atoms(50)

for i = 1:length(unique_perms)
    SRO(unique_perms[i])
end

"For now, it just prints out the SRO parameters for a binary alloy. The ternary
part is already included and fundamentally tested, however the example file is
a binary alloy anyway. As I already said, this file needs further improvements,
it just should serve as a little proof of concept :-)"
println(alpha_ab_array)
println(alpha_aab_array)
println(alpha_abb_array)

"p.s.: prints in Julia are painfully bad"
end
