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
end
