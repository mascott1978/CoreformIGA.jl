module Assembler

import ..Integral
import ..FunctionSpace

function assembleInnerProduct!( i_fc::Integral.FunctionCollectionIntegral,
                                fs_test_fc::FunctionSpace.FunctionCollectionFunctionSpace,
                                fs_trial_fc::FunctionSpace.FunctionCollectionFunctionSpace,
                                disp_strain_mat::Function,
                                test_dof_ids,
                                trial_dof_ids,
                                M )
    for i in 1 : i_fc.quadrature_point_count()
        e_i, xi_i, w_i = i_fc.quadrature_point( i )
        x_i = i_fc.geometry_value( e_i, xi_i )
        weighting_i = i_fc.integral_weight( x_i, e_i, xi_i )

        e_test_i, xi_test_i = fs_test_fc.geometric_map_inversion( x_i, e_i, xi_i )
        M_i = fs_test_fc.global_basis_evaluator( e_test_i, xi_test_i )

        e_trial_i, xi_trial_i = fs_trial_fc.geometric_map_inversion( x_i, e_i, xi_i )
        N_i = fs_trial_fc.global_basis_evaluator( e_trial_i, xi_trial_i )

        for a in 1 : fs_test_fc.global_function_count_on_element( e_test_i )
            a_global_func_id = fs_test_fc.global_function_id_on_element( e_test_i, a )
            a_local_dof_ids = findall( x-> x>0 , test_dof_ids.global_dof_id( a_global_func_id ) )
            temp = cat( test_dof_ids.global_dof_id( a_global_func_id ), dims=1 )
            a_global_dof_ids = temp[ a_local_dof_ids ]
            for b in 1 : fs_trial_fc.global_function_count_on_element( e_trial_i )
                b_global_func_id = fs_trial_fc.global_function_id_on_element( e_trial_i, b )
                b_local_dof_ids = findall( x-> x>0 , trial_dof_ids.global_dof_id( b_global_func_id ) )
                b_local_dofs = cat(trial_dof_ids.global_dof_id( b_global_func_id ), dims = 1 ) 
                b_global_dof_ids = b_local_dofs[ b_local_dof_ids ]
                nodal_M = weighting_i * disp_strain_mat( M_i[ a ] ) * disp_strain_mat( N_i[ b ] )* w_i
                M[ a_global_dof_ids, b_global_dof_ids ] = M[ a_global_dof_ids, b_global_dof_ids ] + nodal_M[ a_local_dof_ids, b_local_dof_ids ]
            end
        end
    end
    return M

end

function assembleProjection!( i_fc::Integral.FunctionCollectionIntegral,
                              fs_test_fc::FunctionSpace.FunctionCollectionFunctionSpace,
                              test_dof_ids, F )

    for i in 1 : i_fc.quadrature_point_count()
        e_i, xi_i, w_i = i_fc.quadrature_point( i )
        x_i = i_fc.geometry_value( e_i, xi_i )
        weighting_i = i_fc.integral_weight( x_i, e_i, xi_i )

        e_test_i, xi_test_i = fs_test_fc.geometric_map_inversion( x_i, e_i, xi_i )
        M_i = fs_test_fc.global_basis_evaluator( e_test_i, xi_test_i )

        for a in 1 : fs_test_fc.global_function_count_on_element( e_test_i )
            a_global_func_id = fs_test_fc.global_function_id_on_element( e_test_i, a )
            a_local_dof_ids = findall( x-> x>0 , test_dof_ids.global_dof_id( a_global_func_id ) )
            temp = cat( test_dof_ids.global_dof_id( a_global_func_id ), dims = 1 )
            a_global_dof_ids =  temp[ a_local_dof_ids ]
            nodal_F = weighting_i * M_i[ a ] * w_i
            F[ a_global_dof_ids ] = F[ a_global_dof_ids ] + nodal_F[ a_local_dof_ids ]
        end
    end
    return F
end

end
