module CutCellDomains1d


struct Layout
    cells::Array{ Array{ Float64, 1 }, 1 }
end

function iterate_domains( layout::Layout, process_domain::Function )
    for e in 1 : size( layout.cells, 1 )
        for cut in 1 : ( size( layout.cells[ e ], 1 ) - 1 )
            x_0 = layout.cells[ e ][ cut ]
            x_1 = layout.cells[ e ][ cut + 1 ]
            process_domain( e, x_0, x_1 )
        end
    end
end

function layout_cut_cell_domains( element_count::Function, d_bc_fc, n_bc_fc, inverse_map )
    cells::Array{ Array{ Float64, 1 }, 1 } = fill([], element_count())
    for e in 1:size( cells, 1 )
        cells[ e ] = [ 0.0, 1.0 ]
    end

    nodes_d = d_bc_fc.nodes( 1 )
    for n in 1 : size( nodes_d, 1 )
        e, x = inverse_map.geometric_map_inversion( nodes_d[ n ], 0, 0.0 )
        if ( abs( x[ 1 ] - 0.0 ) > 1e-10 ) && ( abs( x[ 1 ] - 1.0 ) > 1e-10 )
            append!( cells[ e ], x[ 1 ] )
        end
    end

    nodes_n = n_bc_fc.nodes( 1 )
    for n in 1 : size( nodes_n, 1 )
        e, x = inverse_map.geometric_map_inversion( nodes_n[ n ], 0, 0.0 )
        if ( abs( x[ 1 ] - 0.0 ) > 1e-10 ) && ( abs( x[ 1 ] - 1.0 ) > 1e-10 )
            append!( cells[ e ], x[ 1 ] )
        end
    end

    for e in 1 : size( cells, 1 )
        sort!( cells[ e ] )
    end
    return Layout( cells )
end



end
