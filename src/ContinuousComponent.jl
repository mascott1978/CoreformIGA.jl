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
- `integral`: specify the integrand and quadrature information;
- `geom_field_fc`: offer geometry info. for the inverse mapping and integral weights, etc., in the assembly;
- `test_index_fc`: store the DOFs indexes and related functionalities for the test function space;
- `trial_index_fc`: store the DOFs indexes and related functionalities for the trial function space;
- `test_fs_fc`: test function space for evaluating spline functions and indexing them;
"""
struct FunctionCollection
    integral::Integral.FunctionCollectionIntegral 
    geom_field_fc::Field.FunctionCollection 
    test_index_fc::Index.FunctionCollection
    trial_index_fc::Index.FunctionCollection
    test_fs_fc::FunctionSpace.FunctionCollectionFunctionSpace
    trial_fs_fc::FunctionSpace.FunctionCollectionFunctionSpace
end


function function_collection( ;integral::Integral.FunctionCollectionIntegral, geom_field::Field.FunctionCollection, test_index::Index.FunctionCollection, trial_index::Index.FunctionCollection, test_function_space::FunctionSpace.FunctionCollectionFunctionSpace, trial_function_space::FunctionSpace.FunctionCollectionFunctionSpace )
    return FunctionCollection( integral, geom_field, test_index, trial_index, test_function_space, trial_function_space )
end

end
