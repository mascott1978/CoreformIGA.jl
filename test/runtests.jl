using CoreformIGA
using Test
using LinearAlgebra

@testset "CoreformIGA.jl" begin
    # Write your tests here.
    @testset "BasisBernstein.jl" begin
        @test CoreformIGA.BasisBernstein.B( 1, 0.5 )[ 1 ] == 1.0/2.0
        @test CoreformIGA.BasisBernstein.B( 1, 0.5 )[ 2 ] == 1.0/2.0
        @test CoreformIGA.BasisBernstein.B( 2, 0.5 )[ 1 ] == 0.25
        @test CoreformIGA.BasisBernstein.B( 2, 0.5 )[ 2 ] == 0.5
        @test CoreformIGA.BasisBernstein.B( 2, 0.5 )[ 3 ] == 0.25
        @test CoreformIGA.BasisBernstein.B( 3, 0.5 )[ 1 ] == 1.0/8.0
        @test CoreformIGA.BasisBernstein.B( 3, 0.5 )[ 2 ] == 3.0/8.0
        @test CoreformIGA.BasisBernstein.B( 3, 0.5 )[ 3 ] == 3.0/8.0
        @test CoreformIGA.BasisBernstein.B( 3, 0.5 )[ 4 ] == 1.0/8.0
        @test CoreformIGA.BasisBernstein.dBdxi( 1, 0.5 )[ 1 ] == -1
        @test CoreformIGA.BasisBernstein.dBdxi( 1, 0.5 )[ 2 ] == 1
    end
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
    @testset "BasisSpline.jl" begin
        layout = CoreformIGA.BasisMesh.layout_bspline_1d_uniform_h_max_k( 1, 2, domain = [ 0, 2 ] )
        bm_fc = CoreformIGA.BasisMesh.function_collection( layout )
        bs_fc = CoreformIGA.BasisSpline.function_collection( bm_fc )
        @test bs_fc.global_basis_value( 1, 0.5 ) == [0.5, 0.5]
        @test bs_fc.global_basis_value( 2, 0.4 ) == [0.6, 0.4]
        @test bs_fc.global_basis_parametric_gradient( 1, 0.5 ) == [-1, 1]
        @test bs_fc.global_basis_parametric_gradient( 2, 0.4 ) == [-1, 1]
    end
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
        @test field_fc.field_value( 1, 0.5 ) == [ 1.0, 0.0, 0.0 ]
        @test field_fc.field_value( 2, 0.5 ) == [ 3.0, 0.0, 0.0 ]
        @test field_fc.field_parametric_gradient( 1, 0.5 ) == [ 2.0, 0.0, 0.0 ]
        @test field_fc.field_parametric_gradient( 2, 0.5 ) == [ 2.0, 0.0, 0.0 ]
    end
    @testset "Quadrature.jl" begin
        # qp, qw = CoreformIGA.Quadrature.rule_gauss_legendre_1d( 3 )
        # @test abs( qp[ 1 ] - 0.112702 ) < 1e-5
        # @test abs( qp[ 2 ] - 0.500000 ) < 1e-5
        # @test abs( qp[ 3 ] - 0.887298 ) < 1e-5
        # @test abs( qw[ 1 ] - 0.277778 ) < 1e-5
        # @test abs( qw[ 2 ] - 0.444444 ) < 1e-5
        # @test abs( qw[ 3 ] - 0.277778 ) < 1e-5
        element_count() = 2
        element_degree( e ) = 2
        layout = CoreformIGA.Quadrature.layout_gauss_legendre_1d( element_count, element_degree, #=unused input parameter=# 0 )
        @test abs( layout.quadrature_points[ 1 ] - 0.112702 ) < 1e-5
        @test abs( layout.quadrature_points[ 2 ] - 0.500000 ) < 1e-5
        @test abs( layout.quadrature_points[ 3 ] - 0.887298 ) < 1e-5
        @test abs( layout.quadrature_weights[ 1 ] - 0.277778 ) < 1e-5
        @test abs( layout.quadrature_weights[ 2 ] - 0.444444 ) < 1e-5
        @test abs( layout.quadrature_weights[ 3 ] - 0.277778 ) < 1e-5
        @test abs( layout.quadrature_points[ 4 ] - 0.112702 ) < 1e-5
        @test abs( layout.quadrature_points[ 5 ] - 0.500000 ) < 1e-5
        @test abs( layout.quadrature_points[ 6 ] - 0.887298 ) < 1e-5
        @test abs( layout.quadrature_weights[ 4 ] - 0.277778 ) < 1e-5
        @test abs( layout.quadrature_weights[ 5 ] - 0.444444 ) < 1e-5
        @test abs( layout.quadrature_weights[ 6 ] - 0.277778 ) < 1e-5
        @test layout.element_id_at_quadrature_point[ 1 ] == 1
        @test layout.element_id_at_quadrature_point[ 2 ] == 1
        @test layout.element_id_at_quadrature_point[ 3 ] == 1
        @test layout.element_id_at_quadrature_point[ 4 ] == 2
        @test layout.element_id_at_quadrature_point[ 5 ] == 2
        @test layout.element_id_at_quadrature_point[ 6 ] == 2
        @testset "NonlinearSolver.jl" begin
            # linear
            layout = CoreformIGA.BasisMesh.layout_bspline_1d_uniform_h_max_k( 1, 2, domain = [ 0, 1 ] )
            bm_fc = CoreformIGA.BasisMesh.function_collection( layout )
            bs_fc = CoreformIGA.BasisSpline.function_collection( bm_fc )
            nodes = [ [0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [4.0, 0.0, 0.0] ]
            field_fc = CoreformIGA.Field.function_collection( bm_fc, bs_fc, nodes )
            predictor( x ) = 0.5
            e = 1
            residual( x, xi ) = x - field_fc.field_value( e, xi )
            tangent( xi ) = field_fc.field_parametric_gradient( e, xi )
            solver_linear( K, R ) = ( 1.0 / K ) * R
            norm = LinearAlgebra.norm
            update( xi, delta_xi ) = xi + delta_xi
            solver_nonlinear = CoreformIGA.NonlinearSolver.newtonRaphsonIteration
            x = [0.5, 0, 0]
            xi = solver_nonlinear( x, predictor, residual, tangent, solver_linear, norm, update )
            @test xi == 0.25

            # nonlinear
            layout = CoreformIGA.BasisMesh.layout_bspline_1d_uniform_h_max_k( 2, 2, domain = [ 0, 2 ] )
            bm_fc = CoreformIGA.BasisMesh.function_collection( layout )
            bs_fc = CoreformIGA.BasisSpline.function_collection( bm_fc )
            nodes = [ [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.5, 0.0, 0.0], [2.0, 0.0, 0.0] ]
            field_fc = CoreformIGA.Field.function_collection( bm_fc, bs_fc, nodes )
            predictor( x ) = 0.5
            e = 1
            residual( x, xi ) = x - field_fc.field_value( e, xi )
            tangent( xi ) = field_fc.field_parametric_gradient( e, xi )
            solver_linear( K, R ) = ( 1.0 / K ) * R
            norm = LinearAlgebra.norm
            update( xi, delta_xi ) = xi + delta_xi
            solver_nonlinear = CoreformIGA.NonlinearSolver.newtonRaphsonIteration
            x = [0.5, 0, 0]
            xi = solver_nonlinear( x, predictor, residual, tangent, solver_linear, norm, update )
            @test xi == 0.27924070610240825
        end
        @testset "Geometry.jl" begin
            layout = CoreformIGA.BasisMesh.layout_bspline_1d_uniform_h_max_k( 2, 2, domain = [ 0, 2 ] )
            bm_fc = CoreformIGA.BasisMesh.function_collection( layout )
            bs_fc = CoreformIGA.BasisSpline.function_collection( bm_fc )
            nodes = [ [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.5, 0.0, 0.0], [2.0, 0.0, 0.0] ]
            field_fc = CoreformIGA.Field.function_collection( bm_fc, bs_fc, nodes )
            oneD_inverse_map_fc = CoreformIGA.Geometry.function_collection_map_inversion_1d( bm_fc, field_fc )
            @test oneD_inverse_map_fc.geometric_map_inversion( [0.5, 0, 0 ], 1, 0 ) == ( 1, 0.27924070610240825)
        end
    end
end
