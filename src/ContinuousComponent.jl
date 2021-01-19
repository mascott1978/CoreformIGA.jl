"""
A container to hold a function collection which includes all function spaces and fields used to conduct the integral of boundary conditions or stiffness matrices, etc.
"""
module ContinuousComponent
using ..Index
using ..Field
using ..Quadrature
using ..FunctionSpace
using ..Integral

"""
# Members:
- `integral::Integral`: specify the applied constraint. More specifically, when applying a constraint, we have ``\\lambda * ( u - f )``.
- `index_fc::Index.FunctionCollection`: store the DOFs indexes and related functionalities.
- `geom_field_fc::Field.FunctionCollection`: offer geometry info. for the inverse mapping and integral weights, etc., in the assembly.
- `quadrature_fc::Quadrature.FunctionCollectionQuadrature`: offer quadratures functionalities.
- `fs_N_fc::FunctionSpace.FunctionCollectionFunctionSpace`: offer basis function info
- `fs_dNds_fc::FunctionSpace.FunctionCollectionFunctionSpace`: ( optional ) offer the derivative info of the basis functions.
"""
struct FunctionCollection
    integral::Integral.FunctionCollectionIntegral 
    geom_field_fc::Field.FunctionCollection 
    test_index_fc::Index.FunctionCollection
    trial_index_fc::Index.FunctionCollection
    test_fs_fc::FunctionSpace.FunctionCollectionFunctionSpace
    trial_fs_fc::FunctionSpace.FunctionCollectionFunctionSpace
end


function function_collection( integral::Integral.FunctionCollectionIntegral, geom_field_fc::Field.FunctionCollection, test_index_fc::Index.FunctionCollection, trial_index_fc::Index.FunctionCollection, test_fs_fc::FunctionSpace.FunctionCollectionFunctionSpace, trial_fs_fc::FunctionSpace.FunctionCollectionFunctionSpace )
    return FunctionCollection( integral, geom_field_fc, test_index_fc, trial_index_fc, test_fs_fc, trial_fs_fc )
end

end
