using CoreformIGA
using Test

@testset "Field.jl" begin
    layout = CoreformIGA.BasisMesh.layout_bspline_1d_uniform_h_max_k( 1, 2, domain = [ 0, 2 ] )
    bm_fc = CoreformIGA.BasisMesh.function_collection( layout )
    bs_fc = CoreformIGA.BasisSpline.function_collection( bm_fc )
    nodes = [ [0.0, 0.0, 0.0, 1] [2.0, 0.0, 0.0, 1] [4.0, 0.0, 0.0, 1] ]
    field_fc = CoreformIGA.Field.function_collection( bm_fc, bs_fc, nodes )
    @test field_fc.node_on_element( 1, 1 ) == [0.0, 0.0, 0.0, 1]
    @test field_fc.node_on_element( 1, 2 ) == [2.0, 0.0, 0.0, 1]
    @test field_fc.node_on_element( 2, 1 ) == [2.0, 0.0, 0.0, 1]
    @test field_fc.node_on_element( 2, 2 ) == [4.0, 0.0, 0.0, 1]
    @test field_fc.nodes_on_element( 1 ) == [ [0.0, 0.0, 0.0, 1] [2.0, 0.0, 0.0, 1] ]
    @test field_fc.nodes_on_element( 2 ) == [ [2.0, 0.0, 0.0, 1] [4.0, 0.0, 0.0, 1] ]
    ref_1 = zeros( 3, 1 ) 
    ref_1[ 1 ] = 1
    ref_2 = zeros( 3, 1 ) 
    ref_2[ 1 ] = 3
    ref_3 = zeros( 3, 1 ) 
    ref_3[ 1 ] = 2
    @test field_fc.field_value( 1, [ 0.5 ] ) == ref_1
    @test field_fc.field_value( 2, [ 0.5 ] ) == ref_2
    @test field_fc.field_parametric_gradient( 1, [ 0.5 ] ) == ref_3
    @test field_fc.field_parametric_gradient( 2, [ 0.5 ] ) == ref_3
end
