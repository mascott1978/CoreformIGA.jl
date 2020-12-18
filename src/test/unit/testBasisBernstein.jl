using CoreformIGA
using Test

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
