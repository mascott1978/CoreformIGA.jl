using CoreformIGA
using Test

@testset "Field.jl" begin
    layout = CoreformIGA.BasisMesh.layout_bspline_1d_uniform_h_max_k( 1, 2, domain = [ 0, 2 ] )
    bm_fc = CoreformIGA.BasisMesh.function_collection( layout )
    bs_fc = CoreformIGA.BasisSpline.function_collection( bm_fc )
    nodes = [ [0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [4.0, 0.0, 0.0] ]
    field_fc = CoreformIGA.Field.function_collection( bm_fc, bs_fc, nodes )
    @test field_fc.node_on_element( 1, 1 ) == [0.0, 0.0, 0.0]
    @test field_fc.node_on_element( 1, 2 ) == [2.0, 0.0, 0.0]
    @test field_fc.node_on_element( 2, 1 ) == [2.0, 0.0, 0.0]
    @test field_fc.node_on_element( 2, 2 ) == [4.0, 0.0, 0.0]
    @test field_fc.nodes_on_element( 1 ) == [ [0.0, 0.0, 0.0], [2.0, 0.0, 0.0] ]
    @test field_fc.nodes_on_element( 2 ) == [ [2.0, 0.0, 0.0], [4.0, 0.0, 0.0] ]
    @test field_fc.field_value( 1, [ 0.5 ] ) == [ 1.0, 0.0, 0.0 ]
    @test field_fc.field_value( 2, [ 0.5 ] ) == [ 3.0, 0.0, 0.0 ]
    @test field_fc.field_parametric_gradient( 1, [ 0.5 ] ) == [ 2.0, 0.0, 0.0 ]
    @test field_fc.field_parametric_gradient( 2, [ 0.5 ] ) == [ 2.0, 0.0, 0.0 ]
end
