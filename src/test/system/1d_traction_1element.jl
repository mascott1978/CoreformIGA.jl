using Test
using LinearAlgebra
# stretching of an axial rod, using FRM AL
import CoreformIGA
import BoundaryCondition


@testset "1d_traction_1element.jl" begin
    layout_interior = CoreformIGA.BasisMesh.layout_bspline_1d_uniform_h_max_k( 1, 1 )
    nodes_interior = [0,1]
    quad_rules_interior = [2]
    nodes_constraint_bdry = [0]
    nodes_traction_bdry = [1]
    chi(x) = x >= 0 && x <= 1 ? 1.0 : 1e-12
    penalty_constraint(x) = 1
    constraint(x) = 0
    E(x) = 1
    A(x) = 1
    load(x) = 0
    traction(x) = 1

    dirichlet_bc_layouts 
    dirichlet_bcs_fc = 

    K, M, B, F, G, H = CoreformIGA.FlexRepresentationMethod1d.assemble( layout_interior,
                                                                        nodes_interior,
                                                                        quad_rules_interior,
                                                                        nodes_constraint_bdry,
                                                                        nodes_traction_bdry,
                                                                        chi,
                                                                        penalty_constraint,
                                                                        constraint,
                                                                        E,
                                                                        A,
                                                                        load,
                                                                        traction )

    F_ref = zeros( 2, 1 );
    F_ref[2] = 1.0;
    G_ref = zeros( 1, 1 );
    H_ref = zeros( 2, 1 );
    @test K == [1.0 -1.0; -1.0 1.0]
    @test M == [1.0 0.0; 0.0 0.0]
    @test B == [1.0 0.0]
    @test F == F_ref
    @test G == G_ref
    @test H == H_ref

end
