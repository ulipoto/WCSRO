import YAML
data = YAML.load_file("test.yml")
lattice = data["lattice parameters"]
while length(lattice) < 3
    push!(lattice, lattice[1])
end
lattice = Diagonal(lattice)

rep = data["replication"]
uc = [0 0 0; 0.5 0.5 0; 0.5 0 0.5; 0 0.5 0.5]
uc = uc/maximum(rep)
xtrans = collect(0:1/rep[1]:0.9999)
ytrans = collect(0:1/rep[2]:0.9999)
ztrans = collect(0:1/rep[3]:0.9999)
fcoords = []
for i = 1:length(xtrans), j = 1:length(ytrans), k = 1:length(ztrans)
    for m = 1:size(uc)[1]
        push!(fcoords, [uc[m, 1]+xtrans[i], uc[m, 2]+xtrans[j], uc[m, 3]+xtrans[k]])
    end
end
#println(fcoords)
#println(length(fcoords))

comp = data["composition"]
element_dict = Dict()
comp_vec = []
for (idx, val) in enumerate(comp)
    x = round(val[2]/100*length(fcoords), digits = 0)
    push!(comp_vec, Int(x))
end
if sum(comp_vec) > length(fcoords)
    comp_vec[end] -= 1
elseif sum(comp_vec) < length(fcoords)
    comp_vec[end] += 1
end
println(sum(comp_vec))
for (idx, val) in enumerate(comp)
    element_dict[val[1]] = (idx, comp_vec[idx])
end
element_dict
