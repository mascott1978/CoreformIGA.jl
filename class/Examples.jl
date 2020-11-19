
# stretching of an axial rod, using FRM AL


import CoreformIGA
layout_interior = CoreformIGA.BasisMesh.buildBspline1d( 1, 1 )
nodes_interior = [0,1]
quad_rules_interior = [2]

layout_constraint_bdry = CoreformIGA.BasisMesh.buildBspline0d()
nodes_constraint_bdry = [0]

layout_traction_bdry = CoreformIGA.BasisMesh.buildBspline0d()
nodes_traction_bdry = [1]

chi(x) = x >= 0 && x <= 1 ? 1.0 : 1e-12
penalty_constraint(x) = 1
constraint(x) = 0
E(x) = 1
A(x) = 1
load(x) = 0
traction(x) = 1

quadrature_rule_interior,
extraction_operator_interior,
spline_basis_N_interior,
spline_basis_dNdxi_interior,
id_map_interior = CoreformIGA.FlexRepresentationMethod1d.processBasis1d( layout_interior,
                                                                         quad_rules_interior )

node_map_interior,
geom_x_interior,
geom_dxdxi_interior,
closest_parametric_point_interior = CoreformIGA.FlexRepresentationMethod1d.processGeometry1d( layout_interior, nodes_interior,
                                                                                              spline_basis_N_interior,
                                                                                              spline_basis_dNdxi_interior )

function spline_basis_dNdx_interior( e, xi )
    dNdxi = spline_basis_dNdxi_interior( e, xi )
    dxdxi = geom_dxdxi_interior( e, xi )
    return dNdxi * ( 1.0 / dxdxi )
end

quadrature_rule_constraint_bdry,
extraction_operator_constraint_bdry,
spline_basis_N_constraint_bdry,
spline_basis_dNdxi_constraint_bdry,
id_map_constraint_bdry = CoreformIGA.FlexRepresentationMethod1d.processBasis0d( layout_constraint_bdry )

node_map_constraint_bdry,
geom_x_constraint_bdry,
geom_dxdxi_constraint_bdry,
closest_parametric_point_constraint_bdry = CoreformIGA.FlexRepresentationMethod1d.processGeometry0d( nodes_constraint_bdry )

quadrature_rule_traction_bdry,
extraction_operator_traction_bdry,
spline_basis_N_traction_bdry,
spline_basis_dNdx_traction_bdry,
id_map_traction_bdry = CoreformIGA.FlexRepresentationMethod1d.processBasis0d( layout_traction_bdry )

node_map_traction_bdry,
geom_x_traction_bdry,
geom_dxdxi_traction_bdry,
closest_parametric_point_traction_bdry = CoreformIGA.FlexRepresentationMethod1d.processGeometry0d( nodes_traction_bdry )

weight_K( x ) = chi( x ) * E( x ) * A( x )
weight_M( x ) = penalty_constraint( x )
weight_B( x ) = weight_M( x )
weight_body_force( x ) = load( x ) * chi( x )
weight_traction_force( x ) = traction( x )
weight_G( x ) = constraint( x )
weight_H( x ) = penalty_constraint( x ) * constraint( x )

u_func_n = layout_interior.func_n
lambda_func_n = layout_constraint_bdry.func_n

K = zeros( u_func_n, u_func_n )
M = zeros( u_func_n, u_func_n )
B = zeros( lambda_func_n, u_func_n )
F = zeros( u_func_n, 1 )
G = zeros( lambda_func_n, 1 )
H = zeros( u_func_n, 1 )

K = CoreformIGA.Assemble.assembleInnerProduct!( quadrature_rule_interior,
                                                geom_x_interior,
                                                geom_dxdxi_interior,
                                                weight_K,
                                                closest_parametric_point_interior,
                                                spline_basis_dNdx_interior,
                                                id_map_interior,
                                                closest_parametric_point_interior,
                                                spline_basis_dNdx_interior,
                                                id_map_interior, K )

M = CoreformIGA.Assemble.assembleInnerProduct!( quadrature_rule_constraint_bdry,
                                                geom_x_constraint_bdry,
                                                geom_dxdxi_constraint_bdry,
                                                weight_M,
                                                closest_parametric_point_interior,
                                                spline_basis_N_interior,
                                                id_map_interior,
                                                closest_parametric_point_interior,
                                                spline_basis_N_interior,
                                                id_map_interior, M )

B = CoreformIGA.Assemble.assembleInnerProduct!( quadrature_rule_constraint_bdry,
                                                geom_x_constraint_bdry,
                                                geom_dxdxi_constraint_bdry,
                                                weight_B,
                                                closest_parametric_point_constraint_bdry,
                                                spline_basis_N_constraint_bdry,
                                                id_map_constraint_bdry,
                                                closest_parametric_point_interior,
                                                spline_basis_N_interior,
                                                id_map_interior, B )

F = CoreformIGA.Assemble.assembleProjection!( quadrature_rule_interior,
                                              geom_x_interior,
                                              geom_dxdxi_interior,
                                              weight_body_force,
                                              closest_parametric_point_interior,
                                              spline_basis_N_interior,
                                              id_map_interior, F )

F = CoreformIGA.Assemble.assembleProjection!( quadrature_rule_traction_bdry,
                                              geom_x_traction_bdry,
                                              geom_dxdxi_traction_bdry,
                                              weight_traction_force,
                                              closest_parametric_point_interior,
                                              spline_basis_N_interior,
                                              id_map_interior, F )

G = CoreformIGA.Assemble.assembleProjection!( quadrature_rule_constraint_bdry,
                                              geom_x_constraint_bdry,
                                              geom_dxdxi_constraint_bdry,
                                              weight_G,
                                              closest_parametric_point_constraint_bdry,
                                              spline_basis_N_constraint_bdry,
                                              id_map_constraint_bdry, G )

H = CoreformIGA.Assemble.assembleProjection!( quadrature_rule_constraint_bdry,
                                              geom_x_constraint_bdry,
                                              geom_dxdxi_constraint_bdry,
                                              weight_H,
                                              closest_parametric_point_interior,
                                              spline_basis_N_interior,
                                              id_map_interior, H )
