"""
A container to hold a function collection which includes all function spaces and fields used to conduct the integral of boundary conditions or stiffness matrices, etc.
"""
module ContinuousComponent
using ..Index
using ..Field
using ..Quadrature
using ..FunctionSpace

"""
# Members:
- `f::Function`: specify the applied constraint. More specifically, when applying a constraint, we have ``\\lambda * ( u - f )``.
- `index_fc::Index.FunctionCollection`: store the DOFs indexes and related functionalities.
- `geom_field_fc::Field.FunctionCollection`: offer geometry info. for the inverse mapping and integral weights, etc., in the assembly.
- `quadrature_fc::Quadrature.FunctionCollectionQuadrature`: offer quadratures functionalities.
- `fs_N_fc::FunctionSpace.FunctionCollectionFunctionSpace`: offer basis function info
- `fs_dNds_fc::FunctionSpace.FunctionCollectionFunctionSpace`: ( optional ) offer the derivative info of the basis functions.
"""
struct FunctionCollection
    f::Function 
    index_fc::Index.FunctionCollection
    geom_field_fc::Field.FunctionCollection 
    quadrature_fc::Quadrature.FunctionCollectionQuadrature
    fs_N_fc::FunctionSpace.FunctionCollectionFunctionSpace
    fs_dNds_fc::FunctionSpace.FunctionCollectionFunctionSpace

    function FunctionCollection( f::Function,
                                 index_fc::Index.FunctionCollection,
                                 geom_field_fc::Field.FunctionCollection,
                                 quadrature_fc::Quadrature.FunctionCollectionQuadrature,
                                 fs_N_fc::FunctionSpace.FunctionCollectionFunctionSpace )
        new( f, index_fc, geom_field_fc, quadrature_fc, fs_N_fc )
    end

    function FunctionCollection( f::Function,
                                 index_fc::Index.FunctionCollection,
                                 geom_field_fc::Field.FunctionCollection,
                                 quadrature_fc::Quadrature.FunctionCollectionQuadrature,
                                 fs_N_fc::FunctionSpace.FunctionCollectionFunctionSpace,
                                 fs_dNds_fc::FunctionSpace.FunctionCollectionFunctionSpace)
        new( f, index_fc, geom_field_fc, quadrature_fc, fs_N_fc, fs_dNds_fc )
    end
end


function function_collection( f::Function, index_fc::Index.FunctionCollection, geom_field_fc::Field.FunctionCollection, quadrature_fc::Quadrature.FunctionCollectionQuadrature, fs_N_fc::FunctionSpace.FunctionCollectionFunctionSpace )
    return FunctionCollection( f, index_fc, geom_field_fc, quadrature_fc, fs_N_fc )
end


function function_collection( f::Function, index_fc::Index.FunctionCollection, geom_field_fc::Field.FunctionCollection, quadrature_fc::Quadrature.FunctionCollectionQuadrature, fs_N_fc::FunctionSpace.FunctionCollectionFunctionSpace, fs_dNds_fc::FunctionSpace.FunctionCollectionFunctionSpace )
    return FunctionCollection( f, index_fc, geom_field_fc, quadrature_fc, fs_N_fc, fs_dNds_fc )
end

end
