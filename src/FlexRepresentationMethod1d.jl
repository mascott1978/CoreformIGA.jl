module FlexRepresentationMethod1d

using LinearAlgebra
using Plots
using ..BasisSpline
using ..QuadratureGauss
using ..Field


struct quad_point
    e::Int
    qp::Float64
    qw::Float64
end

function solve( deg, elem_n, quad_rules, flex_domain, p_cad, cad_domain, E, A, load, traction, p_u, constraint )

    layout = BasisSpline.buildUniformHMaxK( deg, elem_n, domain = flex_domain )
    nodes = BasisSpline.nodesEquallySpaced( layout )
    println( "layout: ", layout )
    println( "nodes: ", nodes )

    # callbacks
    function X( e, xi )
        return Field.F( layout, nodes, e, xi )
    end

    function dXdxi( e, xi )
        return Field.dFdxi( layout, nodes, e, xi )
    end

    function N( e, xi )
        return BasisSpline.N( layout, e, xi )
    end

    function dNdxi( e, xi )
        return BasisSpline.dNdxi( layout, e, xi )
    end

    function zero_vector()
        return zeros( layout.func_n, 1 )
    end

    function zero_matrix()
        return zeros( layout.func_n, layout.func_n )
    end

    function id_map( e, a )
        return layout.EG[ e, a ]
    end

    # quadrature points
    qp, qw, IE = QuadratureGauss.legendreLayout( quad_rules )

    quad_points = Array{ quad_point, 1 }( undef, size( qp )[ 1 ] )
    for i in 1:size( qp )[ 1 ]
        quad_points[ i ] = quad_point( IE[ i ], qp[ i ], qw[ i ] )
    end

    # geometry processing
    function closest_point( x )
        return Field.closestPoint( layout, nodes, x )
    end

    function chi( x )
        if x >= cad_domain[ 1 ] && x <= cad_domain[ 2 ]
            return 1.0
        else
            return p_cad
        end
    end


    # assemble linear systems
    K = zero_matrix()
    M = zero_matrix()
    B = zero_vector()
    F = zero_vector()
    G = zeros( 1 )
    H = zero_vector()
    K = assembleK!( quad_points, dNdxi, X, dXdxi, chi, id_map, E, A, K )
    M = assembleM!( closest_point, cad_domain, N, id_map, p_u, M )
    B = assembleB!( closest_point, cad_domain, N, id_map, B )
    F = assembleF!( quad_points, closest_point, cad_domain, N, X, dXdxi, chi, id_map, load, traction, F )
    G = assembleG!( closest_point, cad_domain, X, constraint, G )
    H = assembleH!( closest_point, cad_domain, N, X, id_map, p_u, constraint, H )

    #println( "K: ", K )
    #println( "M: ", M )
    #println( "B: ", B )
    #println( "F: ", F )
    #println( "G: ", G )
    #println( "H: ", H )

    #Uzawa iteration
    count = 30
    d_curr = zero_vector()
    lambda_curr = zeros( 1 )
    for i in 1:count
        K = K + p_u .* M
        rhs_1 = ( F + p_u .* H - lambda_curr .* B  )
        d_curr = K \ rhs_1
        rhs_2 = p_u .* ( transpose(B) * d_curr - G)
        lambda_curr = lambda_curr + rhs_2
        #println( "d_curr: ", d_curr )
        #println( "lambda_curr: ", lambda_curr )
    end

    println( "FINAL d: ", d_curr )
    println( "FINAL lambda: ", lambda_curr )

    plt = Plots.plot()
    for e in 1:elem_n
        Field.plot!( plt, layout, nodes, d_curr, e, Field.F )
    end
    plt
end

function assembleK!( quad_points, dNdxi, X, dXdxi, chi, id_map, E, A, K )
    for i in quad_points
        dNdxi_i = dNdxi( i.e, i.qp )
        X_i = X( i.e, i.qp )
        dXdxi_i = dXdxi( i.e, i.qp )
        p_cad = chi( X_i )
        dNdx_i = dNdxi_i .* ( 1.0 / dXdxi_i )
        for a in 1:size( dNdx_i )[ 1 ]
            for b in 1:size( dNdx_i )[ 1 ]
                res = p_cad * E * A * dNdx_i[ a ] * dNdx_i[ b ] * dXdxi_i * i.qw
                K[ id_map( i.e, a ), id_map( i.e, b ) ] += p_cad * E * A * dNdx_i[ a ] * dNdx_i[ b ] * dXdxi_i * i.qw
            end
        end
    end
    return K
end

function assembleM!( closest_point, cad_domain, N, id_map, p_u, M )
    e_left, xi_left = closest_point( cad_domain[ 1 ] )
    N_left = N( e_left, xi_left )
    for a in 1:size( N_left )[ 1 ]
        for b in 1:size( N_left )[ 1 ]
            M[ id_map( e_left, a ), id_map( e_left, b ) ] += p_u * N_left[ a ] * N_left[ b ]
        end
    end
    return M
end

function assembleB!( closest_point, cad_domain, N, id_map, B )
    e_left, xi_left = closest_point( cad_domain[ 1 ] )
    N_left = N( e_left, xi_left )
    for a in 1:size( N_left )[ 1 ]
        B[ id_map( e_left, a ) ] += N_left[ a ]
    end
    return B
end

function assembleF!( quad_points, closest_point, cad_domain, N, X, dXdxi, chi, id_map, load, traction, F )
    for i in quad_points
        N_i = N( i.e, i.qp )
        X_i = X( i.e, i.qp )
        dXdxi_i = dXdxi( i.e, i.qp )
        p_cad = chi( X_i )
        for a in 1:size( N_i )[ 1 ]
            F[ id_map( i.e, a ) ] += p_cad * N_i[ a ] * load( X_i ) * dXdxi_i * i.qw
        end
    end
    e_right, xi_right = closest_point( cad_domain[ 2 ] )
    N_right = N( e_right, xi_right )
    X_right = X( e_right, xi_right )
    for a in 1:size( N_right )[ 1 ]
        F[ id_map( e_right, a ) ] += N_right[ a ] * traction( X_right )
    end
    return F
end

function assembleG!( closest_point, cad_domain, X, constraint, G )
    e_left, xi_left = closest_point( cad_domain[ 1 ] )
    X_left = X( e_left, xi_left )
    G[ 1 ] = constraint( X_left )
    return G
end

function assembleH!( closest_point, cad_domain, N, X, id_map, p_u, constraint, H )
    e_left, xi_left = closest_point( cad_domain[ 1 ] )
    N_left = N( e_left, xi_left )
    X_left = X( e_left, xi_left )
    for a in 1:size( N_left )[ 1 ]
        H[ id_map( e_left, a ) ] += p_u * N_left[ a ] * constraint( X_left )
    end
    return H
end


end

