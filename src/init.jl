"In this sheet, a given .xyz-file is converted into a matrix. Furthermore, all
necessary arrays for calculating pair- and triplet-SRO parameters are built."

using DelimitedFiles
using LinearAlgebra

function init_matrix()

    "All this parameters are currently necessary for subsequent calculations.
    However, values like limit may shall be input variables later on."
    global limit = 5 #calculating limit
    global sequence = []
    global atom_a = "Bla"
    global atom_b = "Ble"
    global atom_c = "Blu"
    global a_counter = 0
    global b_counter = 0
    global c_counter = 0
    global distances = zeros(1)
    global angles = zeros(1)
    global perimeters = zeros(1)

    "The .xyz-file is being converted here:"
    global total_atoms = parse(Int, readline(open("Al3Ti.xyz", "r")))
    global u = (readdlm("Al3Ti.xyz",'\t', skipstart = 2, skipblanks = true))
    global xyz = zeros(total_atoms, 3)
    global dist_matrix = zeros(total_atoms, total_atoms)
    elements = 1
    for i in 1:total_atoms
        l = split(u[i])
        push!(sequence, l[1])
        xyz[i, 1] = parse(Float64, l[2])
        xyz[i, 2] = parse(Float64, l[3])
        xyz[i, 3] = parse(Float64, l[4])

        "Here, the numbers of atoms A, B and C are being counted"
        atom_a = sequence[1]
        if sequence[i] == atom_a
            a_counter += 1
        elseif sequence[i] != atom_a && elements == 1
            atom_b = sequence[i]
            b_counter += 1
            elements = 2
        elseif elements == 2
            if sequence[i] != atom_a && sequence[i] != atom_b
                atom_c = sequence[i]
                c_counter += 1
            else
                b_counter += 1
            end
        end

        "This determines the number of different lattice constants"
        global ka = 1000
        global kb = 1000
        global kc = 1000
        global edge_lengths = 3
        for u in 1:total_atoms
            if xyz[u, 1] != 0 && xyz[u, 1] < ka
                global ka = xyz[u, 1]
            end
            if xyz[u, 2] != 0 && xyz[u, 2] < kb
                global kb = xyz[u, 2]
            end
            if xyz[u, 3] != 0 && xyz[u, 3] < kc
                global kc = xyz[u, 3]
            end
        end
        if ka == kb && ka == kc
            global edge_lengths = 1
        elseif ka == kb || ka == kc || kb == kc
            global edge_lengths = 2
        end
    end

    "A distance matrix will be created. The value in the matrix at position [i,j]
    represents the norm of the distance between atom i and atom j"
    for i in 1:total_atoms
        a = xyz[i, 1]
        b = xyz[i, 2]
        c = xyz[i, 3]
        for j in i:(size(xyz)[1])
            if round(sqrt((a - xyz[j, 1])^2 + (b - xyz[j, 2])^2 + (c - xyz[j, 3])^2), digits = 6) < limit
                global dist_matrix[i, j] = round(sqrt((a - xyz[j, 1])^2 + (b - xyz[j, 2])^2 + (c - xyz[j, 3])^2), digits = 6)
                #println(dist_matrix[i, j])
            end

    "A list of all different occurring distances is created here:"
            for k in 1:size(distances)[1]
                if round(distances[k], digits = 6) == dist_matrix[i, j]
                    break
                elseif k == size(distances)[1]
                    push!(distances, dist_matrix[i, j])
                end
            end
        end
        sort!(distances)
        popfirst!(distances)
    end

    "This code creates a list of all angles between the atoms i, j and k, where
    i is the atom that sits in the corner"
    global angle_matrix = zeros(total_atoms, total_atoms, total_atoms)
    for i = 1:total_atoms, j = i:total_atoms, k in j:total_atoms
        vector_a = [xyz[j, 1]-xyz[i, 1] xyz[j, 2]-xyz[i, 2] xyz[j, 3]-xyz[i, 3]]
        vector_b = [xyz[k, 1]-xyz[i, 1] xyz[k, 2]-xyz[i, 2] xyz[k, 3]-xyz[i, 3]]
        angle_matrix[i, j, k] = (180/pi)*acos(round((dot(vector_a, vector_b)/(norm(vector_a)*norm(vector_b))), digits = 5))
        "Creates a list of all possible angles:"
        for l in 1:size(angles)[1]
            if angles[l] == angle_matrix[i, j, k]
                break
            elseif l == size(angles)[1] && isnan(angle_matrix[i, j, k]) == false
                push!(angles, angle_matrix[i, j, k])
            end
        end
    end
    sort!(angles)
    popfirst!(angles)

    "A triplet matrix, that is sorted by the perimeter of the resulting triangle"
    global triplet_matrix = zeros(total_atoms, total_atoms, total_atoms)
    for i = 1:total_atoms, j = i:total_atoms, k in j:total_atoms
        if dist_matrix[i, j] != 0 && dist_matrix[i, k] != 0 && dist_matrix[j, k] != 0
            triplet_matrix[i, j, k] = dist_matrix[i, j] + dist_matrix[i, k] + dist_matrix[j, k]
        end

        #A list of all possible perimeters is created here:
        for l in 1:size(perimeters)[1]
            if round(perimeters[l], digits = 6) == round(triplet_matrix[i, j, k], digits = 6)
                break
            elseif l == size(perimeters)[1]
                push!(perimeters, round(triplet_matrix[i, j, k], digits = 6))
            end
        end
    end
    sort!(perimeters)
    popfirst!(perimeters)
end
init_matrix()
