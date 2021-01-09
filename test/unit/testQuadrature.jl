using CoreformIGA
using Test

@testset "Quadrature.jl" begin
    # qp, qw = CoreformIGA.Quadrature.rule_gauss_legendre_1d( 3 )
    # @test abs( qp[ 1 ] - 0.112702 ) < 1e-5
    # @test abs( qp[ 2 ] - 0.500000 ) < 1e-5
    # @test abs( qp[ 3 ] - 0.887298 ) < 1e-5
    # @test abs( qw[ 1 ] - 0.277778 ) < 1e-5
    # @test abs( qw[ 2 ] - 0.444444 ) < 1e-5
    # @test abs( qw[ 3 ] - 0.277778 ) < 1e-5
    element_count() = 3
    element_degree( e ) = 2

    do_nothing() = 0
    inverse( x, e_x, xi_x ) = 3, [ 0.4 ]

    inv = CoreformIGA.Geometry.FunctionCollectionMapInversion( do_nothing, do_nothing, do_nothing, do_nothing, inverse )

    constraint(x) = 0
    traction(x) = 1

    dirichlet_bc_layouts  = CoreformIGA.BoundaryCondition.Layout( 1, [[]], [[[]]], [ CoreformIGA.BasisMesh.layout_bspline_0d() ], [ CoreformIGA.Quadrature.layout_gauss_legendre_0d() ], [ CoreformIGA.Geometry.function_collection_map_inversion_1d ], [ constraint ] )
    dirichlet_bcs_fc = CoreformIGA.BoundaryCondition.function_collection( dirichlet_bc_layouts )

    neumann_bc_layouts  = CoreformIGA.BoundaryCondition.Layout( 1, [[ 2.4 ]], [[[]]], [ CoreformIGA.BasisMesh.layout_bspline_0d() ], [ CoreformIGA.Quadrature.layout_gauss_legendre_0d() ], [ CoreformIGA.Geometry.function_collection_map_inversion_1d ], [ traction ] )
    neumann_bcs_fc = CoreformIGA.BoundaryCondition.function_collection( neumann_bc_layouts )

    layout = CoreformIGA.Quadrature.layout_gauss_legendre_1d( element_count, element_degree, inv, dirichlet_bcs_fc, neumann_bcs_fc, #=unused input parameter=# 0 )
    @test abs( layout.quadrature_points[ 1 ][ 1 ] - 0.112702 ) < 1e-5
    @test abs( layout.quadrature_points[ 2 ][ 1 ] - 0.500000 ) < 1e-5
    @test abs( layout.quadrature_points[ 3 ][ 1 ] - 0.887298 ) < 1e-5
    @test abs( layout.quadrature_weights[ 1 ] - 0.277778 ) < 1e-5
    @test abs( layout.quadrature_weights[ 2 ] - 0.444444 ) < 1e-5
    @test abs( layout.quadrature_weights[ 3 ] - 0.277778 ) < 1e-5
    @test abs( layout.quadrature_points[ 4 ][ 1 ] - 0.112702 ) < 1e-5
    @test abs( layout.quadrature_points[ 5 ][ 1 ] - 0.500000 ) < 1e-5
    @test abs( layout.quadrature_points[ 6 ][ 1 ] - 0.887298 ) < 1e-5
    @test abs( layout.quadrature_weights[ 4 ] - 0.277778 ) < 1e-5
    @test abs( layout.quadrature_weights[ 5 ] - 0.444444 ) < 1e-5
    @test abs( layout.quadrature_weights[ 6 ] - 0.277778 ) < 1e-5
    @test abs( layout.quadrature_points[ 7 ][ 1 ] - ( 0.112702 * 0.4 ) ) < 1e-5
    @test abs( layout.quadrature_points[ 8 ][ 1 ] - ( 0.500000 * 0.4 ) ) < 1e-5
    @test abs( layout.quadrature_points[ 9 ][ 1 ] - ( 0.887298 * 0.4 ) ) < 1e-5
    @test abs( layout.quadrature_weights[ 7 ] - ( 0.277778 * 0.4 ) ) < 1e-5
    @test abs( layout.quadrature_weights[ 8 ] - ( 0.444444 * 0.4 ) ) < 1e-5
    @test abs( layout.quadrature_weights[ 9 ] - ( 0.277778 * 0.4 ) ) < 1e-5
    @test abs( layout.quadrature_points[ 10 ][ 1 ] - ( 0.112702 * 0.6 + 0.4 ) ) < 1e-5
    @test abs( layout.quadrature_points[ 11 ][ 1 ] - ( 0.500000 * 0.6 + 0.4 ) ) < 1e-5
    @test abs( layout.quadrature_points[ 12 ][ 1 ] - ( 0.887298 * 0.6 + 0.4 ) ) < 1e-5
    @test abs( layout.quadrature_weights[ 10 ] - ( 0.277778 * 0.6 ) ) < 1e-5
    @test abs( layout.quadrature_weights[ 11 ] - ( 0.444444 * 0.6 ) ) < 1e-5
    @test abs( layout.quadrature_weights[ 12 ] - ( 0.277778 * 0.6 ) ) < 1e-5
    @test layout.element_id_at_quadrature_point[ 1 ] == 1
    @test layout.element_id_at_quadrature_point[ 2 ] == 1
    @test layout.element_id_at_quadrature_point[ 3 ] == 1
    @test layout.element_id_at_quadrature_point[ 4 ] == 2
    @test layout.element_id_at_quadrature_point[ 5 ] == 2
    @test layout.element_id_at_quadrature_point[ 6 ] == 2
    @test layout.element_id_at_quadrature_point[ 7 ] == 3
    @test layout.element_id_at_quadrature_point[ 8 ] == 3
    @test layout.element_id_at_quadrature_point[ 9 ] == 3
    @test layout.element_id_at_quadrature_point[ 10 ] == 3
    @test layout.element_id_at_quadrature_point[ 11 ] == 3
    @test layout.element_id_at_quadrature_point[ 12 ] == 3
end
