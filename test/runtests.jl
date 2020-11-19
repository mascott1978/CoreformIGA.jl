using CoreformIGA
using Test

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
    @testset "BasisSpline.jl" begin
        layout = CoreformIGA.BasisMesh.buildBspline1d( 1, 2 )
        extraction_operator( e ) = layout.ops[ :, ( e - 1 ) * ( layout.degrees[ e ] + 1 ) + 1 : e * ( layout.degrees[ e ] + 1 ) ]
        @test CoreformIGA.BasisSpline.evaluate( 1, 1, 0.5, extraction_operator, CoreformIGA.BasisBernstein.B )[ 1 ] == 1.0/2.0
    end
    @testset "Field.jl" begin
        layout = CoreformIGA.BasisMesh.buildBspline1d( 1, 2 )
        extraction_operator( e ) = layout.ops[ :, ( e - 1 ) * ( layout.degrees[ e ] + 1 ) + 1 : e * ( layout.degrees[ e ] + 1 ) ]
        basis_evaluator( e, xi ) = CoreformIGA.BasisSpline.evaluate( e, layout.degrees[ e ], xi, extraction_operator, CoreformIGA.BasisBernstein.B )
    end
end
