module Precondition
using LinearAlgebra

function precondition( K, r )
    S = scale( K )
    depend_sets = depend_sets_temp = identify_dependency( S * K * ( S' ), r )
    while ~isempty( depend_sets ) # found linearly independent basis
        grouped_depend_sets = group_dependency( depend_sets_temp )
        # orthonormalize
        row_num, col_num = size( K_star )
        G = Matrix(I, row_num, col_num )
        for depend_fns in grouped_depend_sets
            G[depend_fns, depend_fns] = orthonormalize( K, depend_fns  )
        end
        S = G*D
        depend_sets = identify_dependency( S*K*(S'), r )
        depend_sets_temp = append!( depend_sets_temp, depend_sets )
    end
    return S
end

function scale( K )
    D = Diagonal( 1 ./ sqrt.( diag( K ) ) );
end

function identify_dependency( K_star, r )
    depend_sets = []
    row_num, col_num = size( K_star )
    for irow = 1:row_num
        i_depend_set = [ irow ]
        for icol = irow:col_num
            if K_star( irow, icol ) > r
                push!(i_depend_set, icol)
            end
        end
        if length( i_depend_set ) > 1
            push!( depend_sets, i_depend_set )
        end
    end
end

function group_dependency( depend_sets )
    if length( depend_sets ) == 1
        return depend_sets
    else
        i = 1
        while i > 0
            intersected_set_id = []
            for j = ( i + 1 ):length( depend_sets )
                intersection = intersect( depend_sets[ i ], depend_sets[ j ] )
                if ~isempty( intersection )
                    push!( intersected_set_id, j )
                end
            end
            deleteat!( depend_sets, intersected_set_id )
            if length( depend_sets ) > ( i + 1 )
               i = i + 1
            else
                i = -1 # search is done
            end
        end
        return depend_sets
    end
end

function orthonormalize( K, depend_fns  )
    row_num = length( depend_fns )
    S_sigma = scale( K[ depend_fns, depend_fns ] )
    G = Matrix( I, row_num, row_num )
    for irow = 1:row_num
        # orthogonalize
        G[ irow, 1:( irow - 1 ) ] = -K_star[ depend_fns[ irow ], depend_fns[ 1:( irow - 1 ) ] ]
    end
    K_perb = G * K_star * ( G')
    # normalize
    G = G./diag( K_perb )
    return G
end
