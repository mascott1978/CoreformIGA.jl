module Precondition
using LinearAlgebra

function precondition( K, r )
    S = scale( K )
    depend_sets = depend_sets_temp = identify_dependency( S * K * ( S' ), r )
    while ~isempty( depend_sets ) # found linearly independent basis
        grouped_depend_sets = group_dependency( depend_sets_temp )
        # orthonormalize
        for depend_fns in grouped_depend_sets
            S[depend_fns, depend_fns] = orthonormalize( K[ depend_fns, depend_fns ]  )
        end
        depend_sets = identify_dependency( S*K*(S'), r )
        depend_sets_temp = append!( depend_sets_temp, depend_sets )
    end
    return S
end

function scale( K )
    D = diagm( 1 ./ sqrt.( diag( K ) ) );
end

function identify_dependency( K, r )
    depend_sets = []
    row_num, col_num = size( K )
    for irow = 1:row_num
        i_depend_set = [ irow ]
        for icol = ( irow + 1 ):col_num
            if K[ irow, icol ] > r
                push!(i_depend_set, icol )
            end
        end
        if length( i_depend_set ) > 1
            push!( depend_sets, i_depend_set )
        end
    end
    return depend_sets
end

function group_dependency( depend_sets )
    if length( depend_sets ) == 1
        return depend_sets
    else
        anchor = 2
        set_l = length( depend_sets )
        while anchor <= set_l
            intersected_set_id = []
            for j = anchor:set_l
                intersection = intersect( depend_sets[ anchor-1 ], depend_sets[ j ] )
                if ~isempty( intersection )
                    push!( intersected_set_id, j )
                end
            end
            if isempty( intersected_set_id )
                anchor = anchor + 1
            else
                append!( depend_sets[ anchor - 1 ], vcat( depend_sets[ intersected_set_id ]... ) )
                unique!( depend_sets[ anchor - 1 ] )
                sort!( depend_sets[ anchor - 1 ] )
                deleteat!( depend_sets, intersected_set_id )
                set_l = length( depend_sets )
            end
        end
        return depend_sets
    end
end

function orthonormalize( K  )
    row_num = size( K, 1 )
    S_sigma = scale( K )
    K_star = S_sigma * K * ( S_sigma' )
    G = 1.0*Matrix(I, row_num, row_num )
    for irow = 1:row_num
        # orthogonalize
        G[ irow, 1:( irow - 1 ) ] = - K_star[ irow, 1:( irow - 1 ) ]
    end
    # normalize
    S_sigma = G*S_sigma
    S_sigma = S_sigma./sqrt.(diag( S_sigma * K * ( S_sigma' ) ) )
    return S_sigma
end

end
