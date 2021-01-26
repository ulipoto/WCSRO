using LinearAlgebra
using IterTools

uc = [0 0 0; 0.5 0.5 0; 0.5 0 0.5; 0 0.5 0.5]
rep = [2, 2, 2]
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
println(fcoords)
println(length(fcoords))
