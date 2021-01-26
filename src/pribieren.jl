include("sro_calc.jl")

file = "ab.yml" # Please note, that the file has to be in the same folder.
X = initialize(file)

all_perms = all_permutations(sequence(X))
all_sros = []

for i = 1:size(all_perms)[1]
    push!(all_sros, sro(X, all_perms[i]))
end

length(all_sros)

all_vals = []
push!(all_vals, all_sros[1])

for i = 1:length(all_sros)
    c = 0
    for j = 1:length(all_vals)
        if all_vals[j] == all_sros[i]
            c += 1
            break
        end
    end
    if c == 0
        push!(all_vals, all_sros[i])
    end
end
length(all_vals)
all_vals
