module BasisMesh

import ..BasisBernstein

#Comment with what these members are (types)
struct Layout
    domain
    starts
    lengths
    degrees
    smoothnesses
    ops #CHANGE THIS to handle different degrees TODO
    EG #Element to global function map FIXME for varying degree
    func_n
end

struct FunctionCollection
    element_count::Function
    global_function_count::Function
    global_function_count_on_element::Function
    local_function_count_on_element::Function
    global_function_id_on_element::Function
    global_function_ids_on_element::Function
    parametric_map_value::Function
    parametric_map_gradient::Function
    extraction_operator_on_element::Function
    local_basis_value::Function
    local_basis_parametric_gradient::Function
end

function layout_bspline_0d()
    return bspline_uniform_h_max_k( 0, 1, [ 0, 0 ] )
end

function layout_bspline_1d_uniform_h_max_k( p::Integer, elem_n::Integer; domain = [ 0.0, 1.0 ] )
    elem_size = ( domain[ 2 ] - domain[ 1 ] ) / elem_n
    lengths = [ elem_size for i in 1 : elem_n ]
    starts = zeros( elem_n );
    start = domain[ 1 ]
    for e in 1 : elem_n
        starts[ e ] = start
        start += lengths[ e ]
    end
    degrees = [ p for i in 1 : elem_n ]
    smoothnesses = [ p - 1 for i in 1 : elem_n ]
    smoothnesses[ 1 ] = -1
    append!( smoothnesses, -1 )
    ops = zeros( p + 1, ( p + 1 ) * elem_n ) #[ op[1](3x3), op[2](3x3) ... ]
    EG = zeros( Int32, elem_n, p + 1 ) # Function index map [ [1,2,3], [2,3,4], [3,4,5], [4,5,6] ]
    for e in 1 : elem_n
        ops[ :, ( e - 1 ) * ( p + 1 ) + 1 : e * ( p + 1 ) ] = extraction_operator_on_element_uniform_h_max_k( e, p, elem_n )
        EG[ e, : ] = e : e + p
    end
    func_n = elem_n + p
    return Layout( domain, starts, lengths, degrees, smoothnesses, ops, EG, func_n )
end

function function_collection( layout::Layout )
    return FunctionCollection( element_count( layout ),
                               global_function_count( layout ),
                               global_function_count_on_element( layout ),
                               local_function_count_on_element( layout ),
                               global_function_id_on_element( layout ),
                               global_function_ids_on_element( layout ),
                               parametric_map_value( layout ),
                               parametric_map_gradient( layout ),
                               extraction_operator_on_element( layout ), 
                               local_basis_value(),
                               local_basis_parametric_gradient() )
end

function element_count( layout::Layout )
    return element_count() = length( layout.degrees )
end

function global_function_count( layout::Layout )
    return global_function_count() = layout.func_n
end

function global_function_count_on_element( layout::Layout )
    return global_function_count_on_element( e ) = layout.degrees[ e ] + 1
end

function local_function_count_on_element( layout::Layout )
    return local_function_count_on_element( e ) = layout.degrees[ e ] + 1
end

function global_function_id_on_element( layout::Layout )
    return global_function_id_on_element( e, a ) = layout.EG[ e a ]
end

function global_function_ids_on_element( layout::Layout )
    return global_function_ids_on_element( e ) = [ layout.EG[ e, a ] for a in 1 : layout.degrees[ e ] + 1 ]
end

function parametric_map_value( layout::Layout )
    return parametric_map_value( e, xi ) = layout.starts[ e ] + xi * layout.lengths[ e ]
end

function parametric_map_gradient( layout::Layout )
    return parametric_map_gradient( e, xi ) = layout.lengths[ e ]
end

function extraction_operator_on_element( layout::Layout )
    return extraction_operator_on_element( e ) = layout.ops[ :, ( e - 1 ) * ( layout.degrees[ e ] + 1 ) + 1 : e * ( layout.degrees[ e ] + 1 ) ]
end

function local_basis_value()
    return local_basis_value( p, xi ) = BasisBernstein.B( p, xi )
end

function local_basis_parametric_gradient()
    return local_basis_parametric_gradient( p, xi ) = BasisBernstein.dBdxi( p, xi )
end

function extraction_operator_on_element_uniform_h_max_k( e, p, elem_n )
        if p == 0
            return [ 1.0 ];
        elseif p == 1
            return [ 1 0;
                     0 1 ]
        elseif p == 2
            if elem_n == 1
                return [ 1 0 0;
                         0 1 0;
                         0 0 1 ]
            end
            if e == elem_n
                return [ 1/2 0 0;
                         1/2 1 0;
                         0   0 1 ]
            elseif e == 1
                return [ 1 0 0;
                         0 1 1/2;
                         0 0 1/2 ]
            else
                return [ 1/2 0 0;
                         1/2 1 1/2;
                         0   0 1/2 ]
            end
        elseif p == 3
            if elem_n == 1
                return [ 1 0 0 0;
                         0 1 0 0;
                         0 0 1 0;
                         0 0 0 1 ]
            elseif elem_n == 2
                if e == 1
                    return [ 1 0 0 0;
                             0 1 1/2 1/4;
                             0 0 1/2 1/2;
                             0 0 0 1/4 ]
                else
                    return [ 1/4 0 0 0;
                             1/2 1/2 0 0;
                             1/4 1/2 1 0;
                             0 0 0 1 ]
                end
            elseif elem_n == 3
                if e == 1
                    return [ 1 0 0 0;
                             0 1 1/2 1/4;
                             0 0 1/2 7/12;
                             0 0 0 1/6 ]
                elseif e == 2
                    return [ 1/4 0 0 0;
                             7/12 2/3 1/3 1/6;
                             1/6 1/3 2/3 7/12;
                             0 0 0 1/4 ]
                else
                    return [ 1/6 0 0 0;
                             7/12 1/2 0 0;
                             1/4 1/2 1 0;
                             0 0 0 1 ]
                end
            else
                if e == 1
                    return [ 1 0 0 0;
                             0 1 1/2 1/4;
                             0 0 1/2 7/12;
                             0 0 0 1/6 ]
                elseif e == 2
                    return [ 1/4 0 0 0;
                             7/12 2/3 1/3 1/6;
                             1/6 1/3 2/3 2/3;
                             0 0 0 1/6 ]
                elseif e > 2 && e < elem_n - 1
                    return [ 1/6 0 0 0;
                             2/3 2/3 1/3 1/6;
                             1/6 1/3 2/3 2/3;
                             0 0 0 1/6 ]
                elseif e == elem_n - 1
                    return [ 1/6 0 0 0;
                             2/3 2/3 1/3 1/6;
                             1/6 1/3 2/3 7/12;
                             0 0 0 1/4 ]
                else
                    return [ 1/6 0 0 0;
                             7/12 1/2 0 0;
                             1/4 1/2 1 0;
                             0 0 0 1 ]
                end
            end
        else
            error( "invalid degree specified in extraction operator construction.")
        end
    end
end
