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
    element_degree::Function
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
    ops = [ extraction_operator_on_element_uniform_h_max_k( e, p, elem_n ) for e = 1:elem_n ]
    EG = [ [ i for i = e:e+p ]  for e = 1:elem_n ]
    func_n = elem_n + p
    return Layout( domain, starts, lengths, degrees, smoothnesses, ops, EG, func_n )
end

function layout_bspline_0d()
    return layout_bspline_1d_uniform_h_max_k( 0, 1, domain=[ 0.0, 0.0 ] )
end

function layout_uspline_1d( degrees::Array{<:Integer}, smoothnesses::Array{<:Integer}, lengths::Array{<:Real} )
    domain = [ 0.0, sum( lengths ) ]
    starts = [ sum( [ 0.0, lengths... ][ 1:e ] ) for e in 1:length( lengths ) ]

    local_to_global = global_bernstein_index( degrees )
    N = build_spline_basis( degrees, smoothnesses, lengths )
    func_n = size( N, 1 )

    ops = [ N[ vec( mapslices( row -> any( row[ local_to_global( e, 0 ) : local_to_global( e, degrees[e] ) ] .!= 0 ), N, dims = 2 ) ),
               local_to_global( e, 0 ) : local_to_global( e, degrees[e] ) ]
            for e in 1:length( degrees ) ]

    EG = [ findall( i-> i, vec( mapslices(
        row -> any( row[ local_to_global( e, 0 ) : local_to_global( e, degrees[e] ) ] .!= 0 ),
        N, dims = 2 ) ) )
           for e in 1:length( degrees ) ]

    return Layout( domain, starts, lengths, degrees, smoothnesses, ops, EG, func_n )
end

function function_collection( layout::Layout )
    return FunctionCollection( element_count( layout ),
                               element_degree( layout ),
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

function element_degree( layout::Layout )
    return element_degree( e ) = layout.degrees[ e ]
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
    return global_function_id_on_element( e, a ) = layout.EG[ e ][ a ]
end

function global_function_ids_on_element( layout::Layout )
    return global_function_ids_on_element( e ) = [ layout.EG[ e ][ a ] for a in 1 : layout.degrees[ e ] + 1 ]
end

#NOTE This function will not work for BEXT
function parametric_map_value( layout::Layout )
    return parametric_map_value( e, xi ) = layout.starts[ e ] + xi * layout.lengths[ e ]
end

function parametric_map_gradient( layout::Layout )
    return parametric_map_gradient( e, xi ) = layout.lengths[ e ]
end

function extraction_operator_on_element( layout::Layout )
    return extraction_operator_on_element( e ) = layout.ops[ e ]
end

function local_basis_value()
    return local_basis_value( p, xi ) = BasisBernstein.B( p, xi )
end

function local_basis_parametric_gradient()
    return local_basis_parametric_gradient( p, xi ) = BasisBernstein.dBdxi( p, xi )
end

function global_bernstein_index( degrees::Array{<:Integer} )
    return global_bernstein_index( e, i ) = sum( [ p + 1 for p in degrees[1:e-1] ] ) + i + 1
end

function extraction_operator_on_element_uniform_h_max_k( e, p, elem_n )
    if p == 0
        return reshape( [ 1 ], 1, 1 );
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

########################################################
# U-spline construction algorithm
########################################################

function build_uspline_nullspace_matrix( degrees::Array{<:Integer}, smoothnesses::Array{<:Integer}, lengths )
    local_to_global = global_bernstein_index( degrees )
    num_constraints = sum( [ k + 1 for k in smoothnesses ] )
    num_localfuncs = sum( [ p + 1 for p in degrees ] )
    num_elems = length( degrees )
    S = zeros( num_constraints, num_localfuncs );
    constraint_index = 1;
    for i in 1:length( smoothnesses )
        for k in 0:smoothnesses[i]
            if i > 1
                S[ constraint_index, local_to_global( i-1, 0 ) : local_to_global( i-1, degrees[i-1] ) ] =
                    BasisBernstein.dBdxi( degrees[i-1], k, 1.0 ) / ( lengths[ i - 1 ] )^k
            end

            if i <= num_elems
                S[ constraint_index, local_to_global( i, 0 ) : local_to_global( i, degrees[i] ) ] =
                    -BasisBernstein.dBdxi( degrees[i], k, 0.0 ) / ( lengths[ i ] )^k
            end
            constraint_index += 1
        end
    end
    return S
end

function find_spline_function_supports( degrees::Array{<:Integer}, smoothnesses::Array{<:Integer} )
    local_to_global = global_bernstein_index( degrees )
    corners = [];

    function find_support( seed_elem, seed_localfunc )
        # ----- Helper function -----
        # Determines if spline function should extend into next element
        extendIntoNextCell( curr_elem, curr_localfunc ) = curr_elem <= length( degrees ) &&
            degrees[ curr_elem ] - curr_localfunc <= smoothnesses[ curr_elem + 1 ]

        # ----- Helper function -----
        # Determines how many bernstein coefficients the spline function needs
        # from the next element based on already active smoothness constraints
        nextCellExtension( curr_elem, curr_localfunc ) = smoothnesses[ curr_elem + 1 ] - ( degrees[ curr_elem ] - curr_localfunc )

        # ------ find_support function body -----
        curr_localfunc = seed_localfunc
        curr_elem = seed_elem
        while extendIntoNextCell( curr_elem, curr_localfunc )
            curr_localfunc = nextCellExtension( curr_elem, curr_localfunc )
            curr_elem = curr_elem + 1
            append!( corners, local_to_global( curr_elem, curr_localfunc ) )
        end

        return [ local_to_global( seed_elem, seed_localfunc ) local_to_global( curr_elem, curr_localfunc ) ]
    end

    supports = zeros( Integer, 0, 2 )
    for elem in 1:length( degrees )
        for localfunc in 0 : degrees[elem]
            if local_to_global( elem, localfunc ) in corners
                continue
            end
            append!( corners, local_to_global( elem, localfunc ) )

            supports = vcat( supports, find_support( elem, localfunc ) )
        end
    end
    return supports
end

function compute_spline_function_coefficients( S, support )
    # Remove columns outside support, and any zero rows
    Ssub = S[ vec( mapslices( col -> any( col[ support[1]:support[2] ] .!= 0 ), S, dims = 2 ) ), support[1]:support[2] ]
    # Add the row [ 1 0 0 0 ... ] to the matrix to change nullspace problem into a full-rank system
    Ssub = vcat( Ssub, hcat( [1.0], zeros( 1, length( support[1]:support[2] ) - 1 ) ) )
    # rhs is [ 0 0 ... 0 1 ]'
    b = vcat( zeros( size( Ssub, 1 ) - 1, 1 ), [1.0] )
    # Solve for the coefficients
    coeffs = zeros( size( S, 2 ), 1 )
    coeffs[ support[1]:support[2], 1 ] = Ssub\b
    return coeffs
end

function build_spline_basis( degrees, smoothnesses, lengths )
    S = build_uspline_nullspace_matrix( degrees, smoothnesses, lengths )
    supports = find_spline_function_supports( degrees, smoothnesses )

    N = zeros( size( S, 2 ), size( supports, 1 ) )
    for i in 1:size( supports, 1 )
        N[:,i] = compute_spline_function_coefficients( S, supports[ i, : ] )
    end

    # Normalize N with least-squares
    c = N'*N \ ( N'*ones( size( S, 2 ), 1 ) )
    return c.*N'
end

end
