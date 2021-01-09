using CoreformIGA
using Test
using LinearAlgebra

@testset "NonlinearSolver.jl" begin
    # linear
    layout = CoreformIGA.BasisMesh.layout_bspline_1d_uniform_h_max_k( 1, 2, domain = [ 0, 1 ] )
    bm_fc = CoreformIGA.BasisMesh.function_collection( layout )
    bs_fc = CoreformIGA.BasisSpline.function_collection( bm_fc )
    nodes = [ [0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [4.0, 0.0, 0.0] ]
    field_fc = CoreformIGA.Field.function_collection( bm_fc, bs_fc, nodes )
    predictor( x ) = [ 0.5 ]
    e = 1
    residual( x, xi ) = x - field_fc.field_value( e, xi )
    tangent( xi ) = field_fc.field_parametric_gradient( e, xi )
    solver_linear( K, R ) = [ ( 1.0 / K ) * R ]
    norm = LinearAlgebra.norm
    update( xi, delta_xi ) = xi + delta_xi
    solver_nonlinear = CoreformIGA.NonlinearSolver.newtonRaphsonIteration
    x = [0.5, 0, 0]
    xi = solver_nonlinear( x, predictor, residual, tangent, solver_linear, norm, update )
    @test xi == [ 0.25 ]

    # nonlinear
    layout = CoreformIGA.BasisMesh.layout_bspline_1d_uniform_h_max_k( 2, 2, domain = [ 0, 2 ] )
    bm_fc = CoreformIGA.BasisMesh.function_collection( layout )
    bs_fc = CoreformIGA.BasisSpline.function_collection( bm_fc )
    nodes = [ [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.5, 0.0, 0.0], [2.0, 0.0, 0.0] ]
    field_fc = CoreformIGA.Field.function_collection( bm_fc, bs_fc, nodes )
    predictor( x ) = [ 0.5 ]
    e = 1
    residual( x, xi ) = x - field_fc.field_value( e, xi )
    tangent( xi ) = field_fc.field_parametric_gradient( e, xi )
    solver_linear( K, R ) = [ ( 1.0 / K ) * R ]
    norm = LinearAlgebra.norm
    update( xi, delta_xi ) = xi + delta_xi
    solver_nonlinear = CoreformIGA.NonlinearSolver.newtonRaphsonIteration
    x = [0.5, 0, 0]
    xi = solver_nonlinear( x, predictor, residual, tangent, solver_linear, norm, update )
    @test xi == [ 0.27924070610240825 ]
end
