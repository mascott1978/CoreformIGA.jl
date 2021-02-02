using CoreformIGA
using LinearAlgebra
using Test

@testset "Field.jl" begin
    epsilon = 0.05
    K = 1.0*Matrix( I, 2, 2 )
    K[1, 2] = 1 - epsilon^2
    K[2, 1] = 1 - epsilon^2

    r = 0.9
    depend_sets = CoreformIGA.Precondition.identify_dependency( K, r )
    @test length( depend_sets ) == 1
    @test depend_sets[ 1 ] == [ 1, 2]

    depend_sets = [ [ 1, 2], [ 2, 3 ], [ 4, 5], [ 6, 7], [1, 7] ]
    depend_sets = CoreformIGA.Precondition.group_dependency( depend_sets )
    @test length( depend_sets ) == 2
    @test depend_sets[ 1 ] == [ 1, 2, 3, 6, 7 ]
    @test depend_sets[ 2 ] == [ 4, 5 ]

    K_perp = CoreformIGA.Precondition.orthonormalize( K )
    @test isapprox( K_perp, [1.0 0.0; -14.11560529666 14.15098275354], atol = 1e-10 )

    S = CoreformIGA.Precondition.precondition( K, r )
    @test isapprox( S*K*(S'), [1.0 0.0; 0 1], atol = 1e-10 )

end
