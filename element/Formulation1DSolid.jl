module Formulation1DSolid

struct FunctionCollection
    dis_strain_mat::Function
    dofn::Function
end

function function_collection()
    return FunctionCollection(dis_strain_mat(), dofn())
end

function dis_strain_mat()
    return dis_strain_mat( basis_value ) =  ones( 1, 1 ) * basis_value
end

function dofn()
    return dofn() = 1
end

end
