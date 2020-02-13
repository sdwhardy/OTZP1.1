# Main A* algorithm
function as_Astar(start,goal,domain_data)#node
    as_infCsts(domain_data)# sets all F costs to Infinity
    start=domain_data[start.num]
    goal=domain_data[goal.num]
    openQ=Array{Tuple{Float32,Float32,Int32},1}()
    current=node()
    current=start
    current.closedQ=true
    current.openQ=true
    current.G_cost=0
    current.H_cost=lof_pnt2pnt_dist(start.xy,goal.xy)#calculates straight line distance between points
    current.F_cost=current.G_cost+current.H_cost
    current.parent=current
    current.goal=goal.num
    goal.goal=start.num
    while (current.num != goal.num)
        as_openNodes(current.edges,openQ,goal,domain_data)#adds, updates or removes from open Q
        current=domain_data[openQ[1][3]]
        current.closedQ=true
        deleteat!(openQ,1)
    end
    return current
end

# sets all F costs to Infinity
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

#adds, updates or removes from open Q
function as_openNodes(edgs,oQ,goal,nds)
    for edg in edgs
        if (nds[edg.head].openQ == false && nds[edg.head].closedQ == false)
            as_add2oQ(edg,oQ,nds,goal)
        elseif (nds[edg.head].openQ == true && nds[edg.head].closedQ == false)
            oQ=as_updateOQ(edg,nds,oQ)
        else
            oQ=as_removeOQ(edg,oQ)
        end
    end
    sort!(oQ)
end

#removes node from open Q
function as_removeOQ(edg,oQ)
    for (indx,tpl) in enumerate(oQ)
        if(edg.head == tpl[3])
            deleteat!(oQ,indx)
        end
    end
    return oQ
end

#Updates existing node in Open Q
function as_updateOQ(edg,nds,oQ)
    if (edg.lngth+nds[edg.tail].G_cost < nds[edg.head].G_cost)
        nds[edg.head].parent=nds[edg.tail]
        nds[edg.head].G_cost=edg.lngth+nds[edg.tail].G_cost
        nds[edg.head].F_cost=nds[edg.head].H_cost+nds[edg.head].G_cost
        for (indx,tpl) in enumerate(oQ)
            if (edg.head == tpl[3])
                oQ[indx]=(nds[edg.head].F_cost,nds[edg.head].H_cost,edg.head)
            end
        end
    else
    end
    return oQ
end

#Adds a new node to the Open Q
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
