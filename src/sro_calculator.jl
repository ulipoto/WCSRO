include("init.jl")

"This file uses the constructed matrices from init.jl to calculate the SRO para-
-meters. The shuffling of the given atom sequence shall be done in another file
to preserve some clarity."

function SRO()
    global alpha_ab_array = []
    global alpha_ac_array = []
    global alpha_bc_array = []
    current_sequence = sequence

    #Necessary parameters for pair interactions
    n_aa = 0
    n_ab = 0
    n_ac = 0
    n_bb = 0
    n_bc = 0
    n_cc = 0

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

    "This following code is used for binary alloys:"
    if c_counter == 0
        for m = 1:total_atoms
            n_counter = 0
            for n = m:total_atoms
                if dist_matrix[m, n] != 0 && dist_matrix[m, n] < distances[edge_lengths + 1] && current_sequence[m] != current_sequence[n]
                    n_ab += 1
                    n_counter += 1
                end
            end
            alpha_ab = 1 - n_ab/((a_counter/(a_counter + b_counter))*(b_counter/(a_counter + b_counter))*total_atoms*n_counter)
            push!(alpha_ab_array, alpha_ab)
        end

        #Triplet part
        for i = 1:total_atoms, j = i:total_atoms, k = j:total_atoms
            if triplet_matrix[i, j, k] != 0
                if current_sequence[i] == atom_a
                    if current_sequence[j] == atom_a
                        if current_sequence[k] == atom_a
                            n_aaa += 1
                        else
                            n_aab += 1
                        end
                    else
                        if current_sequence[k] == atom_a
                            n_aab += 1
                        else
                            n_abb += 1
                        end
                    end
                else
                    if current_sequence[j] == atom_a
                        if current_sequence[k] == atom_a
                            n_aab += 1
                        else
                            n_abb += 1
                        end
                    else
                        if current_sequence[k] == atom_a
                            n_abb += 1
                        else
                            n_bbb += 1
                        end
                    end
                end
            end
        end
    else
    "Here, the pair part in case of a ternary alloy begins:"
        for m = 1:total_atoms, n in m:total_atoms
            if dist_matrix[m, n] != 0 && dist_matrix[m, n] < distances[edge_lengths + 1]
                if current_sequence[m] == atom_a
                    if current_sequence[n] == atom_a
                        n_aa += 1
                    elseif current_sequence[n] == atom_b
                        n_ab += 1
                    else
                        n_ac += 1
                    end
                elseif current_sequence[m] == atom_b
                    if current_sequence[n] == atom_a
                        n_ab += 1
                    elseif current_sequence[n] == atom_b
                        n_bb += 1
                    else
                        n_bc += 1
                    end
                elseif current_sequence[m] == atom_c
                    if current_sequence[n] == atom_a
                        n_ac += 1
                    elseif current_sequence[n] == atom_b
                        n_bc += 1
                    else
                        n_cc += 1
                    end
                end
            end
        end
        alpha_ab = 1 - n_ab/((a_counter/(a_counter + b_counter + c_counter))*(b_counter/(a_counter + b_counter + c_counter))*total_atoms*12)
        push!(alpha_ab_array, alpha_ab)
        alpha_ac = 1 - n_ac/((a_counter/(a_counter + b_counter + c_counter))*(c_counter/(a_counter + b_counter + c_counter))*total_atoms*12)
        push!(alpha_ac_array, alpha_ac)
        alpha_bc = 1 - n_bc/((b_counter/(a_counter + b_counter + c_counter))*(c_counter/(a_counter + b_counter + c_counter))*total_atoms*12)
        push!(alpha_bc_array, alpha_bc)
    end
    #Triplet part
    for i = 1:total_atoms, j = i:total_atoms, k = j:total_atoms
    #println(alpha_ab_array)
    #println(alpha_ac_array)
    #println(alpha_bc_array)
end
