module ContinuousComponent
using ..Index
using ..Field
using ..Quadrature
using ..FunctionSpace

struct FunctionCollection
    f::Function #f is a function to specify the values of constraints. More specifically, when applying a constraint, we have \lambda * ( u - f ).
    index_fc::Index.FunctionCollection
    geom_field_fc::Field.FunctionCollection # used for Integral
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
