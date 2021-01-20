using CoreformIGA
using Test
using LinearAlgebra

@testset "Geometry.jl" begin
    layout = CoreformIGA.BasisMesh.layout_bspline_1d_uniform_h_max_k( 2, 2, domain = [ 0, 2 ] )
    bm_fc = CoreformIGA.BasisMesh.function_collection( layout )
    bs_fc = CoreformIGA.BasisSpline.function_collection( bm_fc )
    nodes = [ [0.0, 0.0, 0.0, 1] [1.0, 0.0, 0.0, 1] [1.5, 0.0, 0.0, 1] [2.0, 0.0, 0.0, 1] ]
    field_fc = CoreformIGA.Field.function_collection( bm_fc, bs_fc, nodes )
    oneD_inverse_map_fc = CoreformIGA.Geometry.function_collection_map_inversion_1d( bm_fc, field_fc )
    @test oneD_inverse_map_fc.geometric_map_inversion( [0.5, 0, 0], 1, 0 ) == ( 1, [ 0.27924070610240825 ] )
end
