sing Combinatorics

# Necessary functions are stored in here.

"""
This function determines all different atoms and stores them in a dictionary
while counting the respective number of occurrence.
"""
function element_counter(sequence)
    sequence = sequence
    elements = []
    for i = 1:length(sequence)
        push!(elements, sequence[i])
        if allunique(elements) == false
            pop!(elements)
        end
    end
    el_counter = Dict(val => 0 for (idx, val) in enumerate(elements))
    for i = 1:length(sequence)
        for (key, val) in el_counter
            if key == sequence[i]
                el_counter[key] += 1
            end
        end
    end
    return elements, el_counter
end

"""
The next two functions return the corresponding key in the SRO dictionary for
a given element pair or triplet.
NOTE: While input value D should be a dictionary, "a" should be a vector of
length 2 and 3 respectively.
"""
function get_pair_key(D, a)
    perms = collect(permutations(a, 2))
    for i = 1:length(perms)
        if haskey(D, (perms[i][1], perms[i][2])) == true
            return (perms[i][1], perms[i][2])
        end
    end
end

function get_triplet_key(D, a)
    perms = collect(permutations(a, 3))
    for i = 1:length(perms)
        if haskey(D, (perms[i][1], perms[i][2], perms[i][3])) == true
            return (perms[i][1], perms[i][2], perms[i][3])
        end
    end
end

function pbc(xyz_coordinates)
"""
This function creates a matrix with all shortest distances between atoms i and j
whilst taking periodic boundary conditions into account. This is accomplished by
displacing every lattice into all nearest neighbour cells. Hence, every site
is being displaced 26 times (+ the original position).
"""
    # Parameters
    xyz = xyz_coordinates
    xmax = 0
    ymax = 0
    zmax = 0
    dist_mat = zeros(size(xyz)[1], size(xyz)[1]) # distance matrix

    # This array will be necessary later on when we calculated the angles.
    help_xyz = zeros(size(xyz)[1], size(xyz)[1], 3)

    # Determining the size of the supercell. This is necessary for the displacement.
    for i = 1:size(xyz)[1]
        if xyz[i, 1] > xmax
            xmax = xyz[i, 1]
        end
        if xyz[i, 2] > ymax
            ymax = xyz[i, 2]
        end
        if xyz[i, 3] > zmax
            zmax = xyz[i, 3]
        end
    end

    # This array contains every possible site translation in one "direction".
    delta_matrix = [0 0 0; xmax 0 0; 0 ymax 0; 0 0 zmax; xmax ymax 0;
                    -xmax ymax 0; xmax 0 zmax; -xmax 0 zmax; 0 ymax zmax;
                    0 -ymax zmax; xmax ymax zmax; -xmax ymax zmax;
                    xmax -ymax zmax; xmax ymax -zmax]

    # Here, every translation is being checked in positive and negative direction.
    for i = 1:size(xyz)[1]
        for j = i:size(xyz)[1]
            bestdis = 1e20
            for k = 1:size(delta_matrix)[1]
                # testpos and testneg are the new distance vectors in positive
                # and negative direction of the current displacement.
                testpos = xyz[i, :] - (xyz[j, :] + delta_matrix[k, :])
                testneg = xyz[i, :] - (xyz[j, :] - delta_matrix[k, :])
                if norm(testpos) < bestdis && norm(testpos) != 0
                    bestdis = norm(testpos)
                    help_xyz[i, j, :] = testpos
                    if norm(testneg) < bestdis && norm(testneg) != 0
                        bestdis = norm(testneg)
                        help_xyz[i, j, :] = testneg
                    end
                end
            end
            dist_mat[i, j] = bestdis
        end
    end

    # Now, the angles between atoms i, j and k are being calculated, where i is
    # the atom in the "corner".
    angle_matrix = zeros(size(xyz)[1], size(xyz)[1], size(xyz)[1])
    angles = zeros(1)
    for i = 1:size(xyz)[1], j = i:size(xyz)[1], k in j:size(xyz)[1]
        angle_matrix[i, j, k] = (180/pi)*acos(round((dot(help_xyz[i, j, :], help_xyz[i, k, :])/(norm(help_xyz[i, j, :])*norm(help_xyz[i, k, :]))), digits = 5))
        # Creates a list of all possible angles. Just for checking...
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
    return dist_mat, angle_matrix
end
