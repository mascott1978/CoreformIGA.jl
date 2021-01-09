using CoreformIGA
using Test
using LinearAlgebra

@testset "BEXTImport.jl" begin
    filename = "/../json/quadratic_2_element_bext.json"
    filepath = string(@__DIR__, filename)
    io = open(filepath, "r")
    file = read(io, String)
    layout_interior, nodes = CoreformIGA.ImportBEXT.get1dLayoutFromBEXT( file )
    bm_interior_fc = CoreformIGA.BasisMesh.function_collection( layout_interior )

    @test bm_interior_fc.element_count() == 2
    @test bm_interior_fc.element_degree( 1 ) == [ 2 ]
    @test bm_interior_fc.global_function_count() == 4
    @test bm_interior_fc.global_function_count_on_element( 1 ) == 3
    @test bm_interior_fc.local_function_count_on_element( 1 ) == 3
    @test bm_interior_fc.global_function_id_on_element( 2, 1 ) == 2
    @test bm_interior_fc.global_function_ids_on_element( 2 ) == [ 2, 3, 4 ]
end
