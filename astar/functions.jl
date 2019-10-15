function Astar(start,goal,ocn)#node
    openQ=Array{AstarNode,1}()
    closedQ=Array{AstarNode,1}()
    current=AstarNode()
    nd.node=start
    nd.G_cost=0
    nd.H_cost=lof_pnt2pnt_dist(start.xy,goal.xy)
    nd.F_cost=current.G_cost+current.H_cost

end
