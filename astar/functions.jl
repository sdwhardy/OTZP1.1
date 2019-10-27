
function as_Astar(start,goal,domain_data)#node
    as_infCsts(domain_data)
    openQ=Array{Tuple{Float64,Float64,Int64},1}()
    current=node()
    current=start
    current.closedQ=true
    current.openQ=true
    current.G_cost=0
    current.H_cost=lof_pnt2pnt_dist(start.xy,goal.xy)
    current.F_cost=current.G_cost+current.H_cost
    current.parent=current
    current.goal=goal.num
    goal.goal=start.num
    while(current.num != goal.num)
        as_openNodes(current.edges,openQ,goal,domain_data)
        current=domain_data[openQ[1][3]]
        current.closedQ=true
        deleteat!(openQ,1)
    end
    return current
end

function as_infCsts(domn)
    for node in domn
        node.F_cost=Inf
        node.G_cost=Inf
        node.H_cost=Inf
        node.openQ=false
        node.closedQ=false
        node.goal=0
    end
end

function as_openNodes(edgs,oQ,goal,nds)
    for edg in edgs
        if (nds[edg.head].openQ == false && nds[edg.head].closedQ == false)
            as_add2oQ(edg,oQ,nds,goal)
        elseif (nds[edg.head].openQ == true && nds[edg.head].closedQ == false)
            as_updateOQ(edg,nds)
        else
            as_removeOQ(edg,oQ)
        end
    end
    sort!(oQ)
end

function as_removeOQ(edg,oQ)
    for (indx,tpl) in enumerate(oQ)
        if(edg.head == tpl[3])
            deleteat!(oQ,indx)
        end
    end
end

function as_updateOQ(edg,nds)
    if (edg.lngth+nds[edg.tail].G_cost < nds[edg.head].G_cost)
        nds[edg.head].parent=nds[edg.tail]
        nds[edg.head].G_cost=edg.lngth+nds[edg.tail].G_cost
        nds[edg.head].F_cost=nds[edg.head].H_cost+nds[edg.head].G_cost
    else
    end
end

function as_add2oQ(edg,oQ,nds,goal)
    nds[edg.head].openQ=true
    nds[edg.head].goal=goal.goal
    nds[edg.head].parent=nds[edg.tail]
    nds[edg.head].H_cost=lof_pnt2pnt_dist(nds[edg.head].xy,goal.xy)
    nds[edg.head].G_cost=edg.lngth+nds[edg.tail].G_cost
    nds[edg.head].F_cost=nds[edg.head].H_cost+nds[edg.head].G_cost
    dummy_Tpl=(nds[edg.head].F_cost,nds[edg.head].H_cost,nds[edg.head].num)
    push!(oQ,dummy_Tpl)
end
