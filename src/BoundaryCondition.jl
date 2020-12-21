module BoundaryCondition

struct Layout
    bc_num
    nodes
    basis_meshes
    quadrature_layouts
    geom_inv_maps
    values
end  

struct FunctionCollection
    bc_num::Function
    nodes::Function
    basis_mesh::Function
    quadrature_layout::Function
    geom_inv_map::Function
    value::Function
end

function function_collection( layout::Layout )
    return FunctionCollection( bc_num( layout ),
                               nodes( layout ),
                               basis_mesh( layout ),
                               quadrature_layout( layout ),
                               geom_inv_map( layout ),
                               value( layout )
    )
end

function bc_num( layout::Layout )
    return bc_num() = layout.bc_num
end

function nodes( layout::Layout )
    return nodes( i ) = layout.nodes[ i ]
end

function basis_mesh( layout::Layout )
    return basis_mesh( i ) = layout.basis_meshes[ i ]
end

function quadrature_layout( layout::Layout )
    return quadrature_layout( i ) = layout.quadrature_layouts[ i ]
end

function geom_inv_map( layout::Layout )
    return geom_inv_map( i ) = layout.geom_inv_maps[ i ]
end

function value( layout::Layout )
    return value( i ) = layout.values[ i ]
end

end
