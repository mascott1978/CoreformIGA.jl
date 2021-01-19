module Assembler

using ..Integral
using ..FunctionSpace
using ..Index
using ..ContinuousComponent

function assemble( continuous_components )

    dof_n = 0
    for i = 1:length(  continuous_components )
        dof_n = dof_n > ( continuous_components[ i ].test_index_fc.end_id() ) ? dof_n : continuous_components[i].test_index_fc.end_id()
    end
    K = zeros( dof_n, dof_n )
    F = zeros( dof_n, 1 )


    for i = 1:length(  continuous_components )
        if ( cmp( continuous_components[ i ].integral.destination, "LHS_SYM" ) == 0 ) || ( cmp( continuous_components[ i ].integral.destination, "LHS_NONSYM" )== 0 )
            assembleInnerProduct!( continuous_components[ i ], K )
        elseif cmp( continuous_components[ i ].integral.destination, "RHS" ) == 0
            assembleProjection!( continuous_components[ i ], F )
        else
            error("Destination of the integral is note specified correctly!!!")
        end
    end
    return K, F
end

function assembleInnerProduct!( continuous_comp_fc::ContinuousComponent.FunctionCollection, K )
    for i in 1 : continuous_comp_fc.integral.quadrature_point_count()
        e_i, xi_i, w_i = continuous_comp_fc.integral.quadrature_point( i )
        x_i = continuous_comp_fc.geom_field_fc.field_value( e_i, xi_i )
        e_test_i, xi_test_i = continuous_comp_fc.test_fs_fc.geometric_map_inversion( x_i, e_i, xi_i )
        e_trial_i, xi_trial_i = continuous_comp_fc.trial_fs_fc.geometric_map_inversion( x_i, e_i, xi_i )

        K_mat = continuous_comp_fc.integral.integrand( x_i, xi_test_i, e_test_i, xi_trial_i, e_trial_i, w_i, continuous_comp_fc.geom_field_fc, continuous_comp_fc.test_fs_fc, continuous_comp_fc.trial_fs_fc )

        row_ids = continuous_comp_fc.test_index_fc.global_dof_id( continuous_comp_fc.test_fs_fc.global_function_ids_on_element( e_test_i ) )
        col_ids = continuous_comp_fc.trial_index_fc.global_dof_id( continuous_comp_fc.trial_fs_fc.global_function_ids_on_element( e_trial_i ) )

        row_ids = row_ids[:]
        col_ids = col_ids[:]
        K[ row_ids, col_ids ] += K_mat
        if cmp( continuous_comp_fc.integral.destination, "LHS_NONSYM" ) == 0
            K[ col_ids, row_ids ] += K_mat'
        end
    end
end

function assembleProjection!( continuous_comp_fc::ContinuousComponent.FunctionCollection, F )

    for i in 1 : continuous_comp_fc.integral.quadrature_point_count()
        e_i, xi_i, w_i = continuous_comp_fc.integral.quadrature_point( i )
        x_i = continuous_comp_fc.geom_field_fc.field_value( e_i, xi_i )
        e_test_i, xi_test_i = continuous_comp_fc.test_fs_fc.geometric_map_inversion( x_i, e_i, xi_i )

        F_mat = continuous_comp_fc.integral.integrand( x_i, xi_test_i, e_test_i, w_i, continuous_comp_fc.geom_field_fc, continuous_comp_fc.test_fs_fc )
        row_ids = continuous_comp_fc.test_index_fc.global_dof_id( continuous_comp_fc.test_fs_fc.global_function_ids_on_element( e_test_i ) )

        row_ids = row_ids[:]
        F[ row_ids ] += F_mat
    end
end

end
