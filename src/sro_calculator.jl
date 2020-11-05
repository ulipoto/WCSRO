include("init.jl")

"This file uses the constructed matrices from init.jl to calculate the SRO para-
-meters. The shuffling of the given atom sequence shall be done in another file
to preserve some clarity."
init_matrix()
function SRO()
    seq = sequence
    x_a = a_counter/(a_counter + b_counter + c_counter)
    x_b = b_counter/(a_counter + b_counter + c_counter)
    x_c = c_counter/(a_counter + b_counter + c_counter)
    #Necessary parameters for pair interactions
    n_aa = 0
    n_ab = 0
    n_ac = 0
    n_bb = 0
    n_bc = 0
    n_cc = 0
    global alpha_ab_array = []
    global alpha_ac_array = []
    global alpha_bc_array = []
    global n_counter = 0

    #Necessary parameters for triplet interactions
    n_aaa = 0
    n_aab = 0
    n_abb = 0
    n_bbb = 0
    n_aac = 0
    n_acc = 0
    n_bbc = 0
    n_bcc = 0
    n_ccc = 0
    n_abc = 0
    global alpha_aaa_array = []
    global alpha_aab_array = []
    global alpha_abb_array = []
    global alpha_bbb_array = []
    global alpha_aac_array = []
    global alpha_acc_array = []
    global alpha_bbc_array = []
    global alpha_bcc_array = []
    global alpha_ccc_array = []
    global alpha_abc_array = []

    "This following code is used for binary alloys:"
    if c_counter == 0
        for m = 1:total_atoms, n = m:total_atoms
            if dist_matrix[m, n] != 0 && dist_matrix[m, n] < distances[edge_lengths + 1] #&& seq[m] != seq[n]
                n_ab += 1
                n_counter += 1
            end
        end
        alpha_ab = 1 - n_ab/(x_a*x_b*total_atoms*12)
        if alpha_ab == -Inf
            alpha_ab = -1
        end
        push!(alpha_ab_array, alpha_ab)

        #Triplet part, but with angles
        n_counter = 0
        for i = 1:total_atoms, j = i:total_atoms, k = j:total_atoms
            if triplet_matrix[i, j, k] != 0 && angle_matrix[i, j, k] < angles[trunc(Int, (length(angles)*3/4))] && angle_matrix[i, j, k] > angles[trunc(Int, (length(angles)/4))]
                n_counter += 1
                if seq[i] != seq[j]
                    if seq[i] != seq[k]
                        if seq[i] == atom_a
                            n_abb += 1
                        else
                            n_aab += 1
                        end
                    else
                        if seq[i] == atom_a
                            n_aab += 1
                        else
                            n_abb += 1
                        end
                    end
                else
                    if seq[i] != seq[k]
                        if seq[i] == atom_a
                            n_aab += 1
                        else
                            n_abb += 1
                        end
                    else
                        if seq[i] == atom_a
                            n_aaa += 1
                        else
                            n_bbb += 1
                        end
                    end
                end
            end
        end
        println(n_counter)
        alpha_aab = 1 - n_aab/((x_a*x_a*x_b)*3*n_counter)
        if alpha_aab == -Inf
            alpha_aab = -1
        end
        push!(alpha_aab_array, alpha_aab)
        alpha_abb = 1 - n_abb/((x_a*x_b*x_b)*3*n_counter)
        if alpha_abb == -Inf
            alpha_abb = -1
        end
        push!(alpha_abb_array, alpha_abb)
        println(alpha_ab_array)
        println(alpha_aab_array)
        println(alpha_abb_array)

    else #condition for atom ternary alloy

        #Here, the pair part in case of a ternary alloy begins:
        n_counter = 0
        for m = 1:total_atoms, n in m:total_atoms
            if dist_matrix[m, n] != 0 && dist_matrix[m, n] < distances[edge_lengths + 1]
                n_counter += 1
                if seq[m] == atom_a
                    if seq[n] == atom_a
                        n_aa += 1
                    elseif seq[n] == atom_b
                        n_ab += 1
                    else
                        n_ac += 1
                    end
                elseif seq[m] == atom_b
                    if seq[n] == atom_a
                        n_ab += 1
                    elseif seq[n] == atom_b
                        n_bb += 1
                    else
                        n_bc += 1
                    end
                elseif seq[m] == atom_c
                    if seq[n] == atom_a
                        n_ac += 1
                    elseif seq[n] == atom_b
                        n_bc += 1
                    else
                        n_cc += 1
                    end
                end
            end
        end
        alpha_ab = 1 - n_ab/(x_a*x_b*total_atoms*12)
        if alpha_ab == -Inf
            alpha_ab = -1
        end
        push!(alpha_ab_array, alpha_ab)
        alpha_ac = 1 - n_ac/(x_a*x_c*total_atoms*12)
        if alpha_ac == -Inf
            alpha_ac = -1
        end
        push!(alpha_ac_array, alpha_ac)

        alpha_bc = 1 - n_bc/(x_b*x_c*total_atoms*12)
        if alpha_bc == -Inf
            alpha_bc = -1
        end
        push!(alpha_bc_array, alpha_bc)

        #Now, the part for the triplets:
        n_counter = 0
        for i = 1:total_atoms, j = i:total_atoms, k = j:total_atoms
            if triplet_matrix[i, j, k] != 0 && angle_matrix[i, j, k] < angles[trunc(Int, (length(angles)*3/4))] && angle_matrix[i, j, k] > angles[trunc(Int, (length(angles)/4))]
                n_counter += 1
                if seq[i] != seq[j] && seq[i] != seq[k] && seq[j] != seq[k]
                    n_abc += 1
                elseif seq[i] == seq[j] && seq[i] == seq[k]
                    if seq[i] == atom_a
                        n_aaa += 1
                    elseif seq[i] == atom_b
                        n_bbb += 1
                    else
                        n_ccc += 1
                    end
                elseif seq[i] == seq[j] && seq[i] != seq[k]
                    if seq[i] == atom_a
                        if seq[k] == atom_b
                            n_aab += 1
                        else
                            n_aac += 1
                        end
                    elseif seq[i] == atom_b
                        if seq[k] == atom_a
                            n_abb += 1
                        else
                            n_bbc += 1
                        end
                    else
                        if seq[k] == atom_a
                            n_acc += 1
                        else
                            n_bcc += 1
                        end
                    end
                elseif seq[i] != seq[j] && seq[i] == seq[k]
                    if seq[i] == atom_a
                        if seq[j] == atom_b
                            n_aab += 1
                        else
                            n_aac += 1
                        end
                    elseif seq[i] == atom_b
                        if seq[j] == atom_a
                            n_abb += 1
                        else
                            n_bbc += 1
                        end
                    else
                        if seq[j] == atom_a
                            n_acc += 1
                        else
                            n_bcc += 1
                        end
                    end
                else
                    if seq[i] == atom_a
                        if seq[j] == atom_b
                            n_abb += 1
                        else
                            n_acc += 1
                        end
                    elseif seq[i] == atom_b
                        if seq[j] == atom_a
                            n_aab += 1
                        else
                            n_bcc += 1
                        end
                    else
                        if seq[j] == atom_a
                            n_aac += 1
                        else
                            n_bbc += 1
                        end
                    end
                end
            end
        end
        alpha_aaa = 1 - n_aaa/((x_a*x_a*x_a)*n_counter)
        if alpha_aaa == -Inf
            alpha_aaa = -1
        end
        push!(alpha_aaa_array, alpha_aaa)

        alpha_aab = 1 - n_aab/((x_a*x_a*x_b)*3*n_counter)
        if alpha_aab == -Inf
            alpha_aab = -1
        end
        push!(alpha_aab_array, alpha_aab)

        alpha_abb = 1 - n_abb/((x_a*x_b*x_b)*3*n_counter)
        if alpha_abb == -Inf
            alpha_abb = -1
        end
        push!(alpha_abb_array, alpha_abb)

        alpha_bbb = 1 - n_bbb/((x_b*x_b*x_b)*3*n_counter)
        if alpha_bbb == -Inf
            alpha_bbb = -1
        end
        push!(alpha_bbb_array, alpha_bbb)

        alpha_aac = 1 - n_aac/((x_a*x_a*x_c)*3*n_counter)
        if alpha_aac == -Inf
            alpha_aac = -1
        end
        push!(alpha_aac_array, alpha_aac)

        alpha_acc = 1 - n_acc/((x_a*x_c*x_c)*3*n_counter)
        if alpha_acc == -Inf
            alpha_acc = -1
        end
        push!(alpha_acc_array, alpha_acc)

        alpha_bbc = 1 - n_bbc/((x_b*x_b*x_c)*3*n_counter)
        if alpha_bbc == -Inf
            alpha_bbc = -1
        end
        push!(alpha_bbc_array, alpha_bbc)

        alpha_bcc = 1 - n_bcc/((x_b*x_c*x_c)*3*n_counter)
        if alpha_bcc == -Inf
            alpha_bcc = -1
        end
        push!(alpha_bcc_array, alpha_bcc)

        alpha_ccc = 1 - n_ccc/((x_c*x_c*x_c)*n_counter)
        if alpha_ccc == -Inf
            alpha_ccc = -1
        end
        push!(alpha_ccc_array, alpha_ccc)

        alpha_abc = 1 - (n_abc/((x_a*x_b*x_c)*6*n_counter))
        if alpha_abc == -Inf
            alpha_abc = -1
        end
        push!(alpha_abc_array, alpha_abc)
    end
    println(n_counter, " ", n_aaa, " ", n_aab, " ", n_abb, " ", n_bbb, " ", n_abc)
    println(n_counter)
end
SRO()
println(alpha_abc_array)
println(a_counter)
println(b_counter)
println(c_counter)
