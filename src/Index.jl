"""
Handle DOFs index for each continuous index.
"""
module Index

"""
# Members:
- `start_id::Int64`: the start index of all indexes
- `end_id::Int64`: the end index of all indexes
- `ID`: a one-dimensional array formatted as: [ [], [], ... ]. The ``i``th inter array indicates the index of dofs of ``i``th node.
"""
struct Layout
    start_id::Int64
    end_id::Int64
    ID  
end

"""
# Members:
- `start_id::Function`: get the start index
- `end_id::Function`: get the end index
- `global_dof_id::Function`: get the global DOF index of the ``i``th node.
"""
struct FunctionCollection
    start_id::Function
    end_id::Function
    global_dof_id::Function
end

function function_collection( layout::Layout )
    return FunctionCollection( start_id( layout ),
                               end_id( layout ),
                               global_dof_id( layout ) )
end

function global_dof_id( layout::Layout )
    return global_dof_id( i ) = layout.ID[ :, i ]
end

function start_id( layout::Layout )
    return start_id() = layout.start_id
end

function end_id( layout::Layout )
    return end_id() = layout.end_id
end


# # pass basis information into functions to build this.
# u_func_n = bm_interior_fc.global_function_count()
# node_dofn = formulation.dofn()
# u_func_id2dof_ids = zeros( Int, u_func_n, node_dofn )
# u_dof_n = 0
# for i = 1:u_func_n
#     u_func_id2dof_ids[ i, : ] = [ 1:node_dofn; ] .+ u_dof_n
#     u_dof_n = u_dof_n + node_dofn
# end

# lambda_func_n = 0
# lambda_dof_n = 0
# lambda_func_id2dof_ids = Array{Any}( undef, dirichlet_bcs.bc_num() )
# for i = 1:dirichlet_bcs.bc_num()
#     bm_c_bdry = dirichlet_bcs.basis_mesh( i )
#     bm_c_bdry_fc = BasisMesh.function_collection( bm_c_bdry )
#     i_lambda_func_n = bm_c_bdry_fc.global_function_count()
#     i_dofs = dirichlet_bcs.dofs( i )
#     lambda_func_id2dof_ids[ i ] = zeros( Int, i_lambda_func_n, node_dofn )
#     for j = 1:i_lambda_func_n
#         node_i_constraint_num = length( i_dofs[ j ] )
#         lambda_func_id2dof_ids[ i ][ j, i_dofs[ j ] ] = [ 1:node_i_constraint_num; ] .+ lambda_dof_n
#         lambda_dof_n = lambda_dof_n + node_i_constraint_num
#     end
#     lambda_func_n = lambda_func_n + i_lambda_func_n
# end

end
