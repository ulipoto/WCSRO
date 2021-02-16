import YAML
import LinearAlgebra
import Base.length
import IterTools: groupby
using Printf
using Formatting

struct Structure
    lattice::Array{Float64,2}
    species::Array{String}
    fcoords::Array{Float64,2}
    replicate::Array{Int16}
    composition::Dict

    function Structure(lattice, species, fcoords, repl, comp)
        size(fcoords, 1) != length(species) && throw(ArgumentError("Length of species does not match the number of atomic positions"))
        new(lattice, species, fcoords, repl, comp)
    end
end

volume(s::Structure) = LinearAlgebra.det(s.lattice)
fractional(s::Structure) = s.fcoords
cartesian(s::Structure) = s.fcoords * s.lattice
length(s::Structure) = size(s.fcoords, 1)

slice_structure(s::Structure, indices::Array{Int, 1}) = Structure(s.lattice, s.species[indices], s.fcoords[indices, :], s.replicate, s.composition)
order(s::Structure, by=identity) = slice_structure(s, sortperm(s.species, by=by))

function sublattice(s::Structure, element::String)
    !(element in s.species) && throw(ArgumentError("Atomic species is $(element) part of the structure"))
    slice_structure(s, findall(el -> el == element , s.species))
end

exclude(s::Structure, lattice::String) = slice_structure(s, findall(el -> el != lattice , s.species))

function supercell(s::Structure, sort=false)
    scale = s.replicate
    inverse = collect(map(s -> 1.0/s, scale))
    num_sites = size(s.fcoords, 1)
    scaled_positions = hcat((inverse .* s.fcoords[row, :] for row in 1:num_sites)...)
    scale_matrix = [scale[1] 0 0; 0 scale[2] 0; 0 0 scale[3]]
    fcoords = vcat((
        (scaled_positions +
            hcat(repeat([inverse .* [i, j, k]], num_sites)...)
        )'
        for i=0:scale[1]-1
        for j=0:scale[2]-1
        for k=0:scale[3]-1
    )...)
    species = vcat((s.species for _ in 1:prod(scale))...)
    lattice = scale_matrix * s.lattice
    final = Structure(lattice, species, fcoords, s.replicate, s.composition)
    sort ? order(final) : final
end

function to_poscar(io::IO, s::Structure, direct=true, comment="automatically generated", format=".7f")
    @printf(io, "# %s\n", comment)
    @printf(io, "1.0\n")
    row_format(N) = join(("{$i:$format}" for i in 1:N), " ")
    print_row(row) = printfmtln(io, row_format(length(row)), row...)
    foreach(print_row, eachrow(s.lattice))
    species, amount =  T((g[1], length(g)) for g in groupby(identity, s.species))
    println(io, join(species, " "))
    println(io, join(amount, " "))
    println(io, direct ? "direct" : "cartesian")
    foreach(print_row, eachrow( direct ? fractional(s) : cartesian(s) ))
end

"""
The Struct 'Backpack' is a struct, which shall contain all variables that
have to be calculated just once. It serves to get rid of potentially disturbing
global variables.
"""
mutable struct Backpack
    supercell::Structure
    element_dict::Dict ## Dictionary with the element as a key and a tuple of
                       ## an assigned index and the occurrence of the element.
    mole_fractions::Dict ## Dictionary with all mole fractions (normed to sublattice)
    shell_matrix::Array ## Array with all relationshsips between two atoms
    triplet_shell_matrix::Array ## An array containing all three shell distances
                                ## between atoms i, j and k
    alpha::Array ## Pair SRO array
    beta::Array ## Triplet SRO array
    pair_prefactors::Array ##Array containg pair prefactors
    trip_prefactors::Dict ## Array containing the constant parts of the triplet
                            ##SRO parameter
end
