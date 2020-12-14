module ImportBEXT

import ..BasisMesh
import JSON

function bext( json )
    o = JSON.parse( json )
    #print(o)
    #print("dim ", o["dim"], "\n")
    elements = o["elements"]
    coefficients = o["coefficients"]
    #dense_vector_info = coefficients["dense_vector_info"]
    dense_coefficient_vectors = coefficients["dense_coefficient_vectors"]
    extracted_coeff_vec = [d_coeff_vec["components"] for d_coeff_vec in dense_coefficient_vectors]
    element_blocks = elements["element_blocks"]
    element_block_info = elements[ "element_block_info" ]
    p = []
    for block_info in element_block_info
        append!( p, [ [block_info["deg_s"], block_info["deg_t"], block_info["deg_u"] ] for i=1:block_info["num_elems"] ] )
    end
    element_extraction = [ [extracted_coeff_vec[i+1] for i in block["coeff_vector_ids"]] for block in element_blocks]
    element_nodes = [block["node_ids"] for block in element_blocks]
    #print( "elem_nodes", element_nodes, "\n")
    #print( "element_extraction", element_extraction, "\n")
    nodes = o["nodes"]
    return p, element_extraction, element_nodes, nodes
end

function get1dLayoutFromBEXT( json )
    p, ee, en, ns = bext( json )
    elem_n = size(en)[1]
    degrees = [ x[1] for x in p ] # [ [ 1, -1, -1 ], [ 1, -1, -1 ], [ p-s, p-t, p-u ] ... ]
    smoothnesses = [ ] #Blank for BEXT
    ops = [ vcat(transpose.(vecs)...) for vecs in ee ]
    EG = [x .+ 1 for x in en] #EG(e,i ) EG[e][i]
    func_n = size(ns)[1]
    domain = [ 0.0, elem_n ]
    lengths = [ 1 for i=1:elem_n ]
    starts = [ i - 1.0 for i=1:elem_n ]
    return BasisMesh.Layout( domain, starts, lengths, degrees, smoothnesses, ops, EG, func_n )
end
