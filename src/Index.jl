module Index

struct Layout
    start_id::Int64
    end_id::Int64
    ID # format: [ [], [], ... ], a one-dimensional array. Each inter array indicates the index of dofs of that node 
end

struct FunctionCollection
    start_id::Function
    end_id::Function
    global_dof_id::Function
end

function function_collection( layout::Layout )
    return FunctionCollection( start_id( layout ),
                               end_id( layout ),
                               global_dof_id( layout ) )
end

function global_dof_id( layout::Layout )
    return global_dof_id( i ) = layout.ID[ :, i ]
end

function start_id( layout::Layout )
    return start_id() = layout.start_id
end

function end_id( layout::Layout )
    return end_id() = layout.end_id
end

end
