################################## Testing functions ###########################
#**
function lof_layoutEez_expand_testing(ocn,pcc)
    #define as y or x axis major
    lof_xyMajor(ocn,pcc)
    #set area of owpp
    lof_setAreaOwpp_zeroSize(ocn)
    #set range of MV
    lof_mVrng(ocn)
    #line all boundaries
    lof_bndafy(ocn)
    #add owpp nodes
    lof_owppNodefy(ocn)
    #Add all background nodes
    #lof_nodifySparse_zero_size(ocn)
    lof_nodifySparse(ocn)
    #number all nodes
    lof_numNodes(ocn)
    #add edges for owpp
    lof_owppEdgefy(ocn)
    #add all background edges
    #lof_edgeifySparse_zero_size(ocn)
    lof_edgeifySparse(ocn)
    ocn.buses=vcat(ocn.pccs, ocn.owpps)#collects all buses

    ##########Printing
    for value in ocn.pccs
        print(value.num)
        print(" - ")
        println(value.node.gps)
    end
    for value in ocn.owpps
        print(value.num)
        print(" - ")
        println(value.node.gps)
    end
    println("GPS coordinates projected onto cartesian plane.")
    println("Axis transformed.")
    l_edges=length(ocn.discretedom.edges)
    return ocn,l_edges
end

function lof_setAreaOwpp_zeroSize(ocn)
    for indx=1:length(ocn.owpps)
        #area, height and width
        ocn.owpps[indx].zone.neg_height=ocn.sys.mvCl
        ocn.owpps[indx].zone.pos_height=ocn.sys.mvCl
        ocn.owpps[indx].zone.pos_width=ocn.sys.mvCl
        ocn.owpps[indx].zone.neg_width=ocn.sys.mvCl
    end
end

#=
function lof_nodifySparse_zero_size(ocn)
    konst=10#normal operation
    #place nodes on boundary on no go zones
    for ng in ocn.nogos
        for lns in [ng.wbnd,ng.ebnd]
            for ln in lns
                for y_step = ln.ymn:ocn.sys.prec/konst:ln.ymx
                    dummy_node=node()
                    dummy_node.xy.y=y_step
                    dummy_node.xy.x=y_step*ln.m_findx+ln.b_findx
                    push!(ocn.discretedom.nodes,deepcopy(dummy_node))
                    push!(ng.nodes,ocn.discretedom.nodes[length(ocn.discretedom.nodes)])
                end
                dummy_node=node()
                dummy_node.xy.y=ln.ymx
                dummy_node.xy.x=ln.ymx*ln.m_findx+ln.b_findx
                push!(ocn.discretedom.nodes,deepcopy(dummy_node))
                push!(ng.nodes,ocn.discretedom.nodes[length(ocn.discretedom.nodes)])
            end
        end
        #lns=[ng.sbnd,ng.nbnd][1]
        #ln=lns[1]
        for lns in [ng.sbnd,ng.nbnd]
            for ln in lns
                for x_step = ln.xmn:ocn.sys.prec/konst:ln.xmx
                    dummy_node=node()
                    dummy_node.xy.x=x_step
                    dummy_node.xy.y=x_step*ln.m_findy+ln.b_findy
                    push!(ocn.discretedom.nodes,deepcopy(dummy_node))
                    push!(ng.nodes,ocn.discretedom.nodes[length(ocn.discretedom.nodes)])
                end
                dummy_node=node()
                dummy_node.xy.x=ln.xmx
                dummy_node.xy.y=ln.xmx*ln.m_findy+ln.b_findy
                push!(ocn.discretedom.nodes,deepcopy(dummy_node))
                push!(ng.nodes,ocn.discretedom.nodes[length(ocn.discretedom.nodes)])
            end
        end
    end

    #place nodes on boundary of owpp
    for owpp in ocn.owpps
        for lns in [owpp.zone.wbnd,owpp.zone.ebnd]
            for ln in lns
                for y_step = ln.ymn:ocn.sys.prec/konst:ln.ymx
                    dummy_node=node()
                    dummy_node.xy.y=y_step
                    dummy_node.xy.x=y_step*ln.m_findx+ln.b_findx
                    push!(ocn.discretedom.nodes,deepcopy(dummy_node))
                    push!(owpp.zone.nodes,ocn.discretedom.nodes[length(ocn.discretedom.nodes)])
                end
                dummy_node=node()
                dummy_node.xy.y=ln.ymx
                dummy_node.xy.x=ln.ymx*ln.m_findx+ln.b_findx
                push!(ocn.discretedom.nodes,deepcopy(dummy_node))
                push!(owpp.zone.nodes,ocn.discretedom.nodes[length(ocn.discretedom.nodes)])
            end
        end
        for lns in [owpp.zone.sbnd,owpp.zone.nbnd]
            for ln in lns
                for x_step = ln.xmn:ocn.sys.prec/konst:ln.xmx
                    dummy_node=node()
                    dummy_node.xy.x=x_step
                    dummy_node.xy.y=x_step*ln.m_findy+ln.b_findy
                    push!(ocn.discretedom.nodes,deepcopy(dummy_node))
                    push!(owpp.zone.nodes,ocn.discretedom.nodes[length(ocn.discretedom.nodes)])
                end
                dummy_node=node()
                dummy_node.xy.x=ln.xmx
                dummy_node.xy.y=ln.xmx*ln.m_findy+ln.b_findy
                push!(ocn.discretedom.nodes,deepcopy(dummy_node))
                push!(owpp.zone.nodes,ocn.discretedom.nodes[length(ocn.discretedom.nodes)])
            end
        end
    end

    #place pccs in main node array
    for pcc in ocn.pccs
        push!(ocn.discretedom.nodes,pcc.node)
    end
    unique!(ocn.discretedom.nodes)
end
=#
#=
function lof_edgeifySparse_zero_size(ocn)
    konst=10#zero size
    buildEdge = true
    buildSection = true
    for (indx0,nd_tail) in enumerate(ocn.discretedom.nodes[length(ocn.owpps)+1:length(ocn.discretedom.nodes)])
        for indx1=indx0+1:1:length(ocn.discretedom.nodes)
            nd_head=ocn.discretedom.nodes[indx1]
            if (nd_tail != nd_head)
                verln=lof_lineDirection(nd_tail,nd_head)
                str8line=lof_getStr8line(nd_head,nd_tail)
                if verln == true
                    for y=str8line.ymn:ocn.sys.prec/konst:str8line.ymx
                        x=y*str8line.m_findx+str8line.b_findx
                        buildSection,area=lof_test4pnt(x,y,ocn)
                        if (buildSection == false)
                            buildEdge = false
                            @goto no_str8line
                        end
                    end
                else
                    for x=str8line.xmn:ocn.sys.prec/konst:str8line.xmx
                        y=x*str8line.m_findy+str8line.b_findy
                        buildSection,area=lof_test4pnt(x,y,ocn)
                        if (buildSection == false)
                            buildEdge = false
                            @goto no_str8line
                        end
                    end
                end
            end
            if (buildEdge == true)
                lof_addEdge(nd_tail,nd_head,ocn.discretedom.edges)
            end
            @label no_str8line
            buildEdge = true
        end
    end
end
=#
