module ImportBEXT

import JSON

export bext

#We assume that all information is given in Dense vectors
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
end
