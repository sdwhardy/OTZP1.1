
#=
t=as_Astar(start,goal,domain_data)
ppf_printOcnXY(ocean,[t])
start=ocean.owpps[8].node #513(30.57),520(37.92)
#start=ocean.discretedom.nodes[513]
goal=ocean.discretedom.nodes[75]
domain_data=ocean.discretedom.nodes
=#
function as_Astar(start,goal,domain_data)#node
    as_infCsts(domain_data)
    start=domain_data[start.num]
    goal=domain_data[goal.num]
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
    while (current.num != goal.num)
        as_openNodes(current.edges,openQ,goal,domain_data)
        current=domain_data[openQ[1][3]]
        current.closedQ=true
        deleteat!(openQ,1)
###############################################################
        #75-150-134-447-453-513-588
        #=for i=1:length(openQ)

            if openQ[i][3]==513 println(string(openQ[i][1:2])*" -513- "*string(i))
            elseif openQ[i][3]==453 println(string(openQ[i][1:2])*" -453- "*string(i))
            elseif openQ[i][3]==447 println(string(openQ[i][1:2])*" -447- "*string(i))
            elseif openQ[i][3]==134 println(string(openQ[i][1:2])*" -134- "*string(i))
            elseif openQ[i][3]==150 println(string(openQ[i][1:2])*" -150- "*string(i))

#75-83-114-169-170-142-230-201-261-266-520-588

            elseif openQ[i][3]==520 println(string(openQ[i][1:2])*" -520- "*string(i))
            elseif openQ[i][3]==266 println(string(openQ[i][1:2])*" -266- "*string(i))
            elseif openQ[i][3]==261 println(string(openQ[i][1:2])*" -261- "*string(i))
            elseif openQ[i][3]==201 println(string(openQ[i][1:2])*" -201- "*string(i))
            elseif openQ[i][3]==230 println(string(openQ[i][1:2])*" -230- "*string(i))
            elseif openQ[i][3]==142 println(string(openQ[i][1:2])*" -142- "*string(i))
            elseif openQ[i][3]==170 println(string(openQ[i][1:2])*" -170- "*string(i))
            elseif openQ[i][3]==169 println(string(openQ[i][1:2])*" -169- "*string(i))
            elseif openQ[i][3]==114 println(string(openQ[i][1:2])*" -114- "*string(i))
            elseif openQ[i][3]==83 println(string(openQ[i][1:2])*" -83- "*string(i))
            end

        end
        println(current.num != goal.num)=#
        ####################################################################
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
#=
edgs=current.edges
oQ=openQ
nds=domain_data
nds[447]
=#

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

function as_removeOQ(edg,oQ)
    for (indx,tpl) in enumerate(oQ)
        if(edg.head == tpl[3])
            deleteat!(oQ,indx)
        end
    end
    return oQ
end

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
