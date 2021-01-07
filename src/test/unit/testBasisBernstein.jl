using CoreformIGA
using Test

@testset "BasisBernstein.jl" begin
    @test CoreformIGA.BasisBernstein.B1d( 1, 0.5 )[ 1 ] == 1.0/2.0
    @test CoreformIGA.BasisBernstein.B1d( 1, 0.5 )[ 2 ] == 1.0/2.0
    @test CoreformIGA.BasisBernstein.B1d( 2, 0.5 )[ 1 ] == 0.25
    @test CoreformIGA.BasisBernstein.B1d( 2, 0.5 )[ 2 ] == 0.5
    @test CoreformIGA.BasisBernstein.B1d( 2, 0.5 )[ 3 ] == 0.25
    @test CoreformIGA.BasisBernstein.B1d( 3, 0.5 )[ 1 ] == 1.0/8.0
    @test CoreformIGA.BasisBernstein.B1d( 3, 0.5 )[ 2 ] == 3.0/8.0
    @test CoreformIGA.BasisBernstein.B1d( 3, 0.5 )[ 3 ] == 3.0/8.0
    @test CoreformIGA.BasisBernstein.B1d( 3, 0.5 )[ 4 ] == 1.0/8.0
    @test CoreformIGA.BasisBernstein.dBdxi1d( 1, 0.5 )[ 1 ] == -1
    @test CoreformIGA.BasisBernstein.dBdxi1d( 1, 0.5 )[ 2 ] == 1
    # 2d tests
    @test abs( CoreformIGA.BasisBernstein.B( [ 1, 1 ], [ 0.25, 0.3 ] )[ 1 ] - 0.525 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.B( [ 1, 1 ], [ 0.25, 0.3 ] )[ 2 ] - 0.175 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.B( [ 1, 1 ], [ 0.25, 0.3 ] )[ 3 ] - 0.225 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.B( [ 1, 1 ], [ 0.25, 0.3 ] )[ 4 ] - 0.075 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 1 ], [ 0.25, 0.3 ] )[ 1, 1 ] + 0.7 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 1 ], [ 0.25, 0.3 ] )[ 1, 2 ] + 0.75 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 1 ], [ 0.25, 0.3 ] )[ 2, 1 ] - 0.7 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 1 ], [ 0.25, 0.3 ] )[ 2, 2 ] + 0.25 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 1 ], [ 0.25, 0.3 ] )[ 3, 1 ] + 0.3 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 1 ], [ 0.25, 0.3 ] )[ 3, 2 ] - 0.75 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 1 ], [ 0.25, 0.3 ] )[ 4, 1 ] - 0.3 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 1 ], [ 0.25, 0.3 ] )[ 4, 2 ] - 0.25 ) < 1e-12
    # 3d tests
    @test abs( CoreformIGA.BasisBernstein.B( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 1 ] - 0.238875 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.B( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 2 ] - 0.079625 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.B( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 3 ] - 0.20475 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.B( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 4 ] - 0.06825 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.B( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 5 ] - 0.043875 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.B( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 6 ] - 0.014625 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.B( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 7 ] - 0.128625 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.B( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 8 ] - 0.042875 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.B( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 9 ] - 0.11025 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.B( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 10 ] - 0.03675 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.B( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 11 ] - 0.023625 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.B( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 12 ] - 0.007875 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 1, 1 ] + 0.3185 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 1, 2 ] + 0.6825 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 1, 3 ] + 0.3675 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 2, 1 ] - 0.3185 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 2, 2 ] + 0.2275 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 2, 3 ] + 0.1225 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 3, 1 ] + 0.273 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 3, 2 ] - 0.39 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 3, 3 ] + 0.315 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 4, 1 ] - 0.273 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 4, 2 ] - 0.13 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 4, 3 ] + 0.105 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 5, 1 ] + 0.0585 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 5, 2 ] - 0.2925 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 5, 3 ] + 0.0675 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 6, 1 ] - 0.0585 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 6, 2 ] - 0.0975 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 6, 3 ] + 0.0225 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 7, 1 ] + 0.1715 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 7, 2 ] + 0.3675 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 7, 3 ] - 0.3675 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 8, 1 ] - 0.1715 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 8, 2 ] + 0.1225 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 8, 3 ] - 0.1225 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 9, 1 ] + 0.147 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 9, 2 ] - 0.21 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 9, 3 ] - 0.315 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 10, 1 ] - 0.147 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 10, 2 ] - 0.07 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 10, 3 ] - 0.105 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 11, 1 ] + 0.0315 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 11, 2 ] - 0.1575 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 11, 3 ] - 0.0675 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 12, 1 ] - 0.0315 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 12, 2 ] - 0.0525 ) < 1e-12
    @test abs( CoreformIGA.BasisBernstein.dBdxi( [ 1, 2, 1 ], [ 0.25, 0.3, 0.35 ] )[ 12, 3 ] - 0.0225 ) < 1e-12
end
