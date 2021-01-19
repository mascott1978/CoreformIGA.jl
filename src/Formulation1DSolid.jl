module Formulation1DSolid

struct FunctionCollection
    K_mat::Function
    body_force_vec::Function
    dofn::Function
end

function function_collection()
    return FunctionCollection( K_mat(), body_force_vec(), dofn() ) 
end


function K_mat()

    function K_mat( test_fs_fc, trial_fs_fc, geom_field_fc, xi_test_i, e_test_i, xi_trial_i, e_trial_i )
        M_i = test_fs_fc.global_basis_gradient_evaluator( e_test_i, xi_test_i )
        N_i = trial_fs_fc.global_basis_gradient_evaluator( e_trial_i, xi_trial_i )
        return M_i * ( N_i' ) * ( 1.0 / geom_field_fc.field_parametric_gradient( e_test_i, xi_test_i ) )
    end

    return K_mat
end

function body_force_vec()

    function body_force_vec( test_fs_fc, geom_field_fc, xi_test_i, e_test_i )
        M_i = test_fs_fc.global_basis_evaluator( e_test_i, xi_test_i )
        return M_i * geom_field_fc.field_parametric_gradient( e_test_i, xi_test_i )
    end

    return body_force_vec
end

function dofn()
    return dofn() = 1
end

end
