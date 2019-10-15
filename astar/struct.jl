mutable struct AstarNode
    G_cost::Float64
    H_cost::Float64
    F_cost::Float64
    node::node
    edges::Array{edge}
    openQ::Bool
    closedQ::Bool
    parent::AstarNode
    AstarNode()=(x=new(); x.parent=x)
end
