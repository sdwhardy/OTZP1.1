function Astar(start,goal,ocn)#node
    openQ=Array{AstarNode,1}()
    closedQ=Array{AstarNode,1}()
    current=AstarNode()
    current.node=start
    current.G_cost=0
    current.H_cost=lof_pnt2pnt_dist(start.xy,goal.xy)
    current.F_cost=current.G_cost+current.H_cost
    as_openNode(current,openQ,goal)
end

function as_openNode(crnt,oQ,goal)
    for edg in crnt.edges
        if edg.openQ == false && edg.closedQ == false

    end
end
