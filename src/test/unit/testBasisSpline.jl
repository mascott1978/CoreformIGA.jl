using CoreformIGA
using Test

@testset "BasisSpline.jl" begin
    layout = CoreformIGA.BasisMesh.layout_bspline_1d_uniform_h_max_k( 1, 2, domain = [ 0, 2 ] )
    bm_fc = CoreformIGA.BasisMesh.function_collection( layout )
    bs_fc = CoreformIGA.BasisSpline.function_collection( bm_fc )
    @test bs_fc.global_basis_value( 1, [ 0.5 ] ) == [0.5, 0.5]
    @test bs_fc.global_basis_value( 2, [ 0.4 ] ) == [0.6, 0.4]
    @test bs_fc.global_basis_parametric_gradient( 1, [ 0.5 ] ) == [-1, 1]
    @test bs_fc.global_basis_parametric_gradient( 2, [ 0.4 ] ) == [-1, 1]
end
