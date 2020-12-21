using CoreformIGA
using Test

@testset "BasisMesh.jl" begin
    layout = CoreformIGA.BasisMesh.layout_bspline_1d_uniform_h_max_k( 1, 2 )
    bm_fc = CoreformIGA.BasisMesh.function_collection( layout )
    @test bm_fc.element_count() == 2
    @test bm_fc.global_function_count() == 3
    @test bm_fc.global_function_count_on_element( 1 ) == 2
    @test bm_fc.global_function_count_on_element( 2 ) == 2
    @test bm_fc.local_function_count_on_element( 1 ) == 2
    @test bm_fc.local_function_count_on_element( 2 ) == 2
    @test bm_fc.global_function_ids_on_element( 1 )[ 1 ] == 1
    @test bm_fc.global_function_ids_on_element( 1 )[ 2 ] == 2
    @test bm_fc.global_function_ids_on_element( 2 )[ 1 ] == 2
    @test bm_fc.global_function_ids_on_element( 2 )[ 2 ] == 3
    @test bm_fc.parametric_map_value( 1, 0.5 ) == 0.25
    @test bm_fc.parametric_map_value( 2, 0.5 ) == 0.75
    @test bm_fc.parametric_map_gradient( 1, 0.5 ) == 0.5
    @test bm_fc.parametric_map_gradient( 2, 0.5 ) == 0.5
    @test bm_fc.extraction_operator_on_element( 1 )[ 1, 1 ] == 1.0
    @test bm_fc.extraction_operator_on_element( 1 )[ 1, 2 ] == 0.0
    @test bm_fc.extraction_operator_on_element( 1 )[ 2, 1 ] == 0.0
    @test bm_fc.extraction_operator_on_element( 1 )[ 2, 2 ] == 1.0
    @test bm_fc.extraction_operator_on_element( 2 )[ 1, 1 ] == 1.0
    @test bm_fc.extraction_operator_on_element( 2 )[ 1, 2 ] == 0.0
    @test bm_fc.extraction_operator_on_element( 2 )[ 2, 1 ] == 0.0
    @test bm_fc.extraction_operator_on_element( 2 )[ 2, 2 ] == 1.0
    @test isapprox( bm_fc.local_basis_value( 1, 0.5 )[ 1 ], 0.5 )
    @test isapprox( bm_fc.local_basis_value( 1, 0.5 )[ 2 ], 0.5 )
    @test bm_fc.local_basis_parametric_gradient( 1, 0.5 )[ 1 ] == -1.0
    @test bm_fc.local_basis_parametric_gradient( 1, 0.5 )[ 2 ] == 1.0
end

@testset "BasisMesh.jl U-splines" begin
    # Matching a B-spline layout
    bspline_layout = CoreformIGA.BasisMesh.layout_bspline_1d_uniform_h_max_k( 2, 3, domain = [ 0.0, 3.0 ] )
    bspline_fc = CoreformIGA.BasisMesh.function_collection( bspline_layout )
    uspline_layout = CoreformIGA.BasisMesh.layout_uspline_1d( [ 2, 2, 2 ], [ -1, 1, 1, -1 ], [ 1.0, 1.0, 1.0 ] )
    uspline_fc = CoreformIGA.BasisMesh.function_collection( uspline_layout )
    @test bspline_fc.element_count() == uspline_fc.element_count()
    @test bspline_fc.global_function_count() == uspline_fc.global_function_count()
    for i in 1:bspline_fc.element_count()
        @test bspline_fc.element_degree( i ) == uspline_fc.element_degree( i )
        @test bspline_fc.global_function_count_on_element( i ) == uspline_fc.global_function_count_on_element( i )
        @test bspline_fc.local_function_count_on_element( i ) == uspline_fc.local_function_count_on_element( i )
        for f in 1:bspline_fc.global_function_count_on_element( i )
            @test bspline_fc.global_function_ids_on_element( i )[ f ] == uspline_fc.global_function_ids_on_element( i )[ f ]
        end
        @test all( x -> isapprox( x... ), zip( bspline_fc.extraction_operator_on_element( i ), uspline_fc.extraction_operator_on_element( i ) ) )
        @test isapprox( bspline_fc.parametric_map_value( i, 0.5 ), uspline_fc.parametric_map_value( i, 0.5 ) )
        @test isapprox( bspline_fc.parametric_map_gradient( i, 0.5 ), uspline_fc.parametric_map_gradient( i, 0.5 ) )
    end

    # Something with non-unit lengths
    layout = CoreformIGA.BasisMesh.layout_uspline_1d( [ 2, 2 ], [ -1, 1, -1 ], [ 0.5, 1.0 ] )
    fc = CoreformIGA.BasisMesh.function_collection( layout )
    @test fc.global_function_count() == 4
    ops1 = [ 1.0 0.0 0.0
             0.0 1.0 2.0/3.0
             0.0 0.0 1.0/3.0 ]
    ops2 = [ 2.0/3.0 0.0 0.0
             1.0/3.0 1.0 0.0
             0.0 0.0 1.0 ]
    @test all( x -> isapprox( x... ), zip( fc.extraction_operator_on_element( 1 ), ops1 ) )
    @test all( x -> isapprox( x... ), zip( fc.extraction_operator_on_element( 2 ), ops2 ) )
    @test fc.global_function_ids_on_element( 1 )[ 1 ] == 1
    @test fc.global_function_ids_on_element( 1 )[ 2 ] == 2
    @test fc.global_function_ids_on_element( 1 )[ 3 ] == 3
    @test fc.global_function_ids_on_element( 2 )[ 1 ] == 2
    @test fc.global_function_ids_on_element( 2 )[ 2 ] == 3
    @test fc.global_function_ids_on_element( 2 )[ 3 ] == 4


    # One more, slightly more complicated.  Ops to compare against are taken from the igx codebase.
    layout = CoreformIGA.BasisMesh.layout_uspline_1d( [ 3, 3, 2 ], [ -1, 2, 1, -1 ], [ 1.0, 3.0, 2.0 ] )
    fc = CoreformIGA.BasisMesh.function_collection( layout )
    ops1 = [ 1      0      0      0
             0      1   0.75 0.5625
             0      0   0.25  0.375
             0      0      0 0.0625 ]
    ops2 = [ 0.5625      0      0      0
             0.375   0.75      0      0
             0.0625   0.25      1    0.5
             0      0      0    0.5 ]
    ops3 = [ 0.5   0   0
             0.5   1   0
             0   0   1 ]
    @test all( x -> isapprox( x... ), zip( fc.extraction_operator_on_element( 1 ), ops1 ) )
    @test all( x -> isapprox( x... ), zip( fc.extraction_operator_on_element( 2 ), ops2 ) )
    @test all( x -> isapprox( x... ), zip( fc.extraction_operator_on_element( 3 ), ops3 ) )
    @test fc.global_function_ids_on_element( 1 )[ 1 ] == 1
    @test fc.global_function_ids_on_element( 1 )[ 2 ] == 2
    @test fc.global_function_ids_on_element( 1 )[ 3 ] == 3
    @test fc.global_function_ids_on_element( 1 )[ 4 ] == 4
    @test fc.global_function_ids_on_element( 2 )[ 1 ] == 2
    @test fc.global_function_ids_on_element( 2 )[ 2 ] == 3
    @test fc.global_function_ids_on_element( 2 )[ 3 ] == 4
    @test fc.global_function_ids_on_element( 2 )[ 4 ] == 5
    @test fc.global_function_ids_on_element( 3 )[ 1 ] == 4
    @test fc.global_function_ids_on_element( 3 )[ 2 ] == 5
    @test fc.global_function_ids_on_element( 3 )[ 3 ] == 6
end
