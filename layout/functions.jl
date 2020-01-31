function lof_layoutEez_basis()
    #reading data from excel input files
    ocean=eez()
    ocean.id_count=0
    ocean.pccs,ocean.id_count=lof_getPccData(ocean.id_count)
    ocean.owpps,ocean.id_count=lof_getOwppData(ocean.id_count)
    ocean.bndryPnts=lof_getBndData()
    ocean.sys=lof_getSysData()
    ocean.finance=cstD_cfs()
    lof_nogoZones(ocean)

    #transformation from gps to cartesian
    ocean.base=lof_bseCrd(ocean)#find base coordinates
    print("Base coordinates are: ")
    println(ocean.base)
    lof_gps2cartesian4bus(ocean.owpps,ocean.base)#projects owpps onto cartesian plane
    lof_gps2cartesian4bus(ocean.pccs,ocean.base)#projects pccs onto cartesian plane
    lof_gps2cartesian4nodes(ocean.bndryPnts,ocean.base)#projects boundary points onto cartesian plane
    for i=1:ocean.sys.nogoNum
        lof_gps2cartesian4nodes(ocean.nogos[i].bndryPnts,ocean.base)#projects boundary points onto cartesian plane
    end

    lof_transformAxis(ocean)
    return ocean
end

function lof_layoutEez_expand(ocn,pcc)
    #define as y or x axis major
    lof_xyMajor(ocn,pcc)
    #set area of owpp
    lof_setAreaOwpp(ocn)
    #set range of MV
    lof_mVrng(ocn)
    #line all boundaries
    lof_bndafy(ocn)
    #Add all background nodes
    lof_nodifySparse(ocn)
    #number all nodes
    lof_numNodes(ocn)

    #add all background edges
    lof_edgeifySparse(ocn)
    #add edges for owpp
    lof_owppEdgefy(ocn)

    ocn.buses=vcat(ocn.pccs, ocn.owpps)#collects all buses
    #calculate constraint ellipse approximations for owpps
    opt_owppConstraints(ocn)
    #calculate constraint ellipse approximations for nogos
    opt_nogoConstraints(ocn)


    ##########Printing
    for value in ocn.pccs
        print(value.num)
        print(" - ")
        println(value.node.gps)
    end
    for value in ocean.owpps
        print(value.num)
        print(" - ")
        println(value.node.gps)
    end
    println("GPS coordinates projected onto cartesian plane.")
    println("Axis transformed.")
    return ocn
end

function lof_layoutEez_expand_testing(ocn,pcc)
    #define as y or x axis major
    lof_xyMajor(ocn,pcc)
    #set area of owpp
    lof_setAreaOwpp(ocn)
    #set range of MV
    lof_mVrng(ocn)
    #line all boundaries
    lof_bndafy(ocn)
    #Add all background nodes
    lof_nodifySparse(ocn)
    #number all nodes
    lof_numNodes(ocn)

    #add all background edges
    #lof_edgeifySparse(ocn)
    #add edges for owpp
    #lof_owppEdgefy(ocn)

    ocn.buses=vcat(ocn.pccs, ocn.owpps)#collects all buses

    ##########Printing
    for value in ocn.pccs
        print(value.num)
        print(" - ")
        println(value.node.gps)
    end
    for value in ocean.owpps
        print(value.num)
        print(" - ")
        println(value.node.gps)
    end
    println("GPS coordinates projected onto cartesian plane.")
    println("Axis transformed.")
    return ocn
end

function lof_xyMajor(ocn,pcc)
    x_min=pcc.node.xy.x
    x_max=pcc.node.xy.x
    y_min=pcc.node.xy.y
    y_max=pcc.node.xy.y
    for owp in ocn.owpps
        if (owp.node.xy.x>x_max)
            x_max=deepcopy(owp.node.xy.x)
        elseif (owp.node.xy.x<x_min)
            x_min=deepcopy(owp.node.xy.x)
        end
        if (owp.node.xy.y>y_max)
            y_max=deepcopy(owp.node.xy.y)
        elseif (owp.node.xy.y<y_min)
            y_min=deepcopy(owp.node.xy.y)
        end
    end
    if (abs(x_max-x_min)>abs(y_max-y_min))
        ocn.yaxisMajor=false
    else
        ocn.yaxisMajor=true
    end
end

function lof_order2Pcc(ocn,pcc)
    lngth_owpps=Array{Tuple{Float64,Int64},1}()
    ordrd_owpps=Array{bus,1}()
    for (indx,owp) in enumerate(ocn.owpps)
        owp.num=deepcopy(indx)
        push!(lngth_owpps,(deepcopy(lof_pnt2pnt_dist(owp.node.xy,pcc.node.xy)),deepcopy(owp.num)))
    end
    sort!(lngth_owpps, by = x -> x[1])
    for indx in lngth_owpps
        for owp in ocn.owpps
            if (indx[2] == owp.num)
                push!(ordrd_owpps,owp)
            end
        end
    end
    for (i,owp) in enumerate(ordrd_owpps)
        owp.num=deepcopy(i)
        owp.id=owp.num+length(ocn.pccs)-1
    end
    return ordrd_owpps
end

function lof_owppEdgefy(ocn)
    for (ind0,owpp) in enumerate(ocn.owpps)
        owpp.node.num=length(ocn.discretedom.nodes)+1
        push!(ocn.discretedom.nodes,owpp.node)
        for (ind1,periNd) in enumerate(owpp.zone.nodes)
            dummy_edge=edge()
            dummy_edge.tail=owpp.node.num
            dummy_edge.head=periNd.num
            dummy_edge.lngth=0.5#length to connect to any point on perimeter of windfarm
            push!(owpp.node.edges,dummy_edge)
        end
    end
end

function lof_edgeifySparse(ocn)
    konst=3
    buildEdge = true
    buildSection = true
    for (indx0,nd_tail) in enumerate(ocn.discretedom.nodes)
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

function lof_edgeifyOss(oss,ocn,nodes,edges)
    buildEdge = true
    buildSection = true
    nd_tail=oss
    for indx1=1:1:length(nodes)-1
        nd_head=nodes[indx1]
        if (nd_tail != nd_head)
            verln=lof_lineDirection(nd_tail,nd_head)
            str8line=lof_getStr8line(nd_head,nd_tail)
            if verln == true
                for y=str8line.ymn:ocn.sys.prec:str8line.ymx
                    x=y*str8line.m_findx+str8line.b_findx
                    buildSection,area=lof_test4pnt(x,y,ocn)
                    if (buildSection == false)
                        buildEdge = false
                        @goto no_str8line
                    end
                end
            else
                #=
                x=str8line.xmn+ocn.sys.prec
                =#
                for x=str8line.xmn:ocn.sys.prec:str8line.xmx
                    y=x*str8line.m_findy+str8line.b_findy
                    buildSection,area=lof_test4pnt(x,y,ocn)
                    if (buildSection == false)
                        #println("Post x: "*string(buildSection)*" - y= "*string(y)*" - x= "*string(x))
                        buildEdge = false
                        @goto no_str8line
                    end
                end
            end
        end
        if (buildEdge == true)
            #println("Building Edge!!")
            lof_addEdge(nd_tail,nd_head,edges)
        end
        @label no_str8line
        buildEdge = true
    end
end

function lof_pntWithinNG(x,y,area2tst)
    inside=true
    vertices=xy[]
    for pnt in area2tst.bndryPnts
        push!(vertices,deepcopy(pnt.xy))
        vertices[length(vertices)].x=(vertices[length(vertices)].x-x)
        vertices[length(vertices)].y=(vertices[length(vertices)].y-y)
    end
    #println(vertices)
    vertices=vertices[end:-1:1]
    a_i=Float64
    pos=false
    neg=false
    nll=false
    for i=1:1:length(vertices)
        if (i==length(vertices))
            j=1
        else
            j=i+1
        end
        a_i=(vertices[j].x*vertices[i].y)-(vertices[i].x*vertices[j].y)
        if (a_i>0)
            pos=true
        elseif (a_i<0)
            neg=true
        else
            nll=true
        end
    end
    if (pos==true && neg==true) || (nll==true)
        inside=false
    end
    return inside
end
#Complete these 2 functions then test in edges are correct
function lof_pntWithin(x,y,area2tst)
    #return true if inside area
    inside=false
    east_in=false
    east_out=false
    west_in=false
    west_out=false
    north_in=false
    north_out=false
    south_in=false
    south_out=false

    #check north boundary
    #ln=area2tst.nbnd[1]
    if length(area2tst.nbnd) != 0
        for ln in area2tst.nbnd
            #if (ln.xmn<=x && ln.xmx>=x && ln.ymx>=y)
            if (ln.xmn<=x && ln.xmx>=x)
                if ((ln.m_findy*x+ln.b_findy)>y)
                    north_in=true
                elseif ((ln.m_findy*x+ln.b_findy)<=y)
                    north_out=true
                end
            end
        end
    end

    #check south boundary
    #ln=area2tst.sbnd[1]
    if length(area2tst.sbnd) != 0
        for ln in area2tst.sbnd
            #if (ln.xmn<=x && ln.xmx>=x && ln.ymn<=y)
            if (ln.xmn<=x && ln.xmx>=x)
                if ((ln.m_findy*x+ln.b_findy)<y)
                    south_in=true
                elseif ((ln.m_findy*x+ln.b_findy)>=y)
                    south_out=true
                end
            end
        end
    end

    #check east boundary
    #ln=area2tst.ebnd[1]
    if length(area2tst.ebnd) != 0
        for ln in area2tst.ebnd
            #if (ln.ymn<=y && ln.ymx>=y && ln.xmx>=x)
            if (ln.ymn<=y && ln.ymx>=y)
                if ((ln.m_findx*y+ln.b_findx)>x)
                    east_in=true
                elseif ((ln.m_findx*y+ln.b_findx)<=x)
                    east_out=true
                end
            end
        end
    end

    #check west boundary
    #ln=area2tst.wbnd[1]
    if length(area2tst.wbnd) != 0
        for ln in area2tst.wbnd
            #if (ln.ymn<=y && ln.ymx>=y && ln.xmn<=x)
            if (ln.ymn<=y && ln.ymx>=y)
                if ((ln.m_findx*y+ln.b_findx)<x)
                    west_in=true
                elseif ((ln.m_findx*y+ln.b_findx)>=x)
                    west_out=true
                end
            end
        end
    end
    if (((west_in && east_in) || (west_in && north_in) || (west_in && south_in) || (east_in && north_in) || (east_in && south_in) || (north_in && south_in)) && (north_out != true && south_out != true && west_out != true && east_out != true))
        inside=true
    end

    return inside
end

function lof_addEdge(nd_tail,nd_head,Alledges)
    #build and store both directions of edge
    lnth=lof_pnt2pnt_dist(nd_tail.xy,nd_head.xy)
    dummy_edgeAB=edge()
    dummy_edgeBA=edge()

    dummy_edgeAB.tail=nd_tail.num
    dummy_edgeAB.head=nd_head.num
    dummy_edgeAB.lngth=lnth

    dummy_edgeBA.tail=nd_head.num
    dummy_edgeBA.head=nd_tail.num
    dummy_edgeBA.lngth=lnth

    push!(nd_tail.edges,deepcopy(dummy_edgeAB))
    push!(nd_head.edges,deepcopy(dummy_edgeBA))
    push!(Alledges,nd_head.edges[length(nd_head.edges)])
    push!(Alledges,nd_tail.edges[length(nd_tail.edges)])
end
#stage_go
#owpp=ocn.owpps[3]
function lof_test4pnt(x,y,ocn)
    within = false
    outside = true
    area=nogo()
    #check owpps
    for owpp in ocn.owpps
        within=lof_pntWithin(x,y,owpp.zone)
        if (within == true)
            area.ebnd=owpp.zone.ebnd
            area.wbnd=owpp.zone.wbnd
            area.nbnd=owpp.zone.nbnd
            area.sbnd=owpp.zone.sbnd
            area.nodes=owpp.zone.nodes
            outside=false
            @goto pnt_inside
        end
    end

    #check nogo areas
    #ngo = ocn.nogos[2]
    for ngo in ocn.nogos
        within=lof_pntWithinNG(x,y,ngo)
        if (within == true)
            area.ebnd=ngo.ebnd
            area.wbnd=ngo.wbnd
            area.nbnd=ngo.nbnd
            area.sbnd=ngo.sbnd
            area.nodes=ngo.nodes
            outside=false
            @goto pnt_inside
        end
    end
    @label pnt_inside
    return outside,area
end

function lof_lineDirection(nd_tail,nd_head)
    if (abs(nd_tail.xy.y-nd_head.xy.y) <  abs(nd_tail.xy.x-nd_head.xy.x))
        vertLn=false
    else
        vertLn=true
    end
    return vertLn
end


function lof_getStr8line(nd_hd,nd_tl)
    dummy_line=line()
    dummy_line.xmx=max(nd_hd.xy.x,nd_tl.xy.x)
    dummy_line.xmn=min(nd_hd.xy.x,nd_tl.xy.x)
    dummy_line.ymx=max(nd_hd.xy.y,nd_tl.xy.y)
    dummy_line.ymn=min(nd_hd.xy.y,nd_tl.xy.y)

    if (nd_hd.xy.x==nd_tl.xy.x)
        dummy_line.b_findy,dummy_line.m_findy=reverse([[nd_hd.xy.x,nd_tl.xy.x+(10^-5)] ones(2)]\[nd_hd.xy.y,nd_tl.xy.y])
    else
        dummy_line.b_findy,dummy_line.m_findy=reverse([[nd_hd.xy.x,nd_tl.xy.x] ones(2)]\[nd_hd.xy.y,nd_tl.xy.y])
    end

    if (nd_hd.xy.y==nd_tl.xy.y)
        dummy_line.b_findx,dummy_line.m_findx=reverse([[nd_hd.xy.y,nd_tl.xy.y+(10^-5)] ones(2)]\[nd_hd.xy.x,nd_tl.xy.x])
    else
        dummy_line.b_findx,dummy_line.m_findx=reverse([[nd_hd.xy.y,nd_tl.xy.y] ones(2)]\[nd_hd.xy.x,nd_tl.xy.x])
    end
    return dummy_line
end

function lof_posPccs(ocn)
    for pcc in ocn.pccs
        bsf_distance=Inf
        bsf_node=Int32
        for (indx,node) in enumerate(ocn.discretedom.nodes)
            distance=lof_pnt2pnt_dist(pcc.node.xy,node.xy)
            if  distance < bsf_distance
                bsf_distance=deepcopy(distance)
                bsf_node=deepcopy(indx)
            end
        end
        println("bsf: "*string(bsf_node)*"at "*string(bsf_distance))
        ocn.discretedom.nodes[bsf_node]=pcc.node
    end
end

function lof_numNodes(ocn)
    for (indx,nd) in enumerate(ocn.discretedom.nodes)
        nd.num=deepcopy(indx)
    end
end
#ng=ocn.nogos[1]
function lof_nodifySparse(ocn)
    #place nodes on boundary on no go zones
    for ng in ocn.nogos
        for lns in [ng.wbnd,ng.ebnd]
            for ln in lns
                for y_step = ln.ymn:ocn.sys.prec:ln.ymx
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
                for x_step = ln.xmn:ocn.sys.prec:ln.xmx
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
                for y_step = ln.ymn:ocn.sys.prec:ln.ymx
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
                for x_step = ln.xmn:ocn.sys.prec:ln.xmx
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

#=function lof_edgeify(ocn)
    hyp=2*ocn.sys.prec
    for (indx0,nd0) in enumerate(ocn.discretedom.nodes)
        for indx1 = indx0+1:1:length(ocn.discretedom.nodes)
            lnth=lof_pnt2pnt_dist(nd0.xy,ocn.discretedom.nodes[indx1].xy)
            if (lnth < hyp)
                dummy_edgeAB=edge()
                dummy_edgeBA=edge()

                dummy_edgeAB.tail=nd0.num
                dummy_edgeAB.head=ocn.discretedom.nodes[indx1].num
                dummy_edgeAB.lngth=lnth

                dummy_edgeBA.tail=ocn.discretedom.nodes[indx1].num
                dummy_edgeBA.head=nd0.num
                dummy_edgeBA.lngth=lnth

                push!(nd0.edges,deepcopy(dummy_edgeAB))
                push!(ocn.discretedom.nodes[indx1].edges,deepcopy(dummy_edgeBA))
                push!(ocn.discretedom.edges,nd0.edges[length(nd0.edges)])
                push!(ocn.discretedom.edges,ocn.discretedom.nodes[indx1].edges[length(ocn.discretedom.nodes[indx1].edges)])
            end
        end
    end
    for (ind0,owpp) in enumerate(ocn.owpps)
        owpp.node.num=length(ocn.discretedom.nodes)+1
        push!(ocn.discretedom.nodes,owpp.node)
        for (ind1,periNd) in enumerate(owpp.zone.pnts)
            dummy_edge=edge()
            dummy_edge.tail=owpp.node.num
            dummy_edge.head=periNd.num
            dummy_edge.lngth=0.5#length to connect to any point on perimeter of windfarm
            push!(owpp.node.edges,dummy_edge)
        end
    end
end=#

function lof_nogoZones(ocn)
    for i=1:ocn.sys.nogoNum
        nogoT=nogo()
        nogoT.bndryPnts=lof_getNoGoData(i)
        push!(ocn.nogos,deepcopy(nogoT))
    end
end

function lof_mVrng(ocn)
    for owpp in ocn.owpps
        mvrng=cstF_mVrng(5,owpp.mva,owpp.wnd,ocn.finance,ocn.sys.prec)
        owpp.mv_zone.pos_height=deepcopy(mvrng)+owpp.zone.pos_height
        owpp.mv_zone.pos_width=deepcopy(mvrng)+owpp.zone.pos_width
        owpp.mv_zone.neg_height=deepcopy(mvrng)+owpp.zone.neg_height
        owpp.mv_zone.neg_width=deepcopy(mvrng)+owpp.zone.neg_width
        owpp.mv_zone.area=(owpp.mv_zone.pos_height+owpp.mv_zone.neg_height)*(owpp.mv_zone.pos_width+owpp.mv_zone.neg_width)
    end
end

function lof_setAreaOwpp(ocn)
    #find dominant axis
    if (abs(ocn.owpps[1].node.xy.y-ocn.owpps[length(ocn.owpps)].node.xy.y)>abs(ocn.owpps[1].node.xy.x-ocn.owpps[length(ocn.owpps)].node.xy.x))
        y_axisDominant=true
    else
        y_axisDominant=false
    end

    if y_axisDominant==true
        #set areas, heights and widths
        dif=abs(ocn.owpps[1].node.xy.y-ocn.owpps[2].node.xy.y)
        ocn.owpps[1].zone.area=ocn.owpps[1].mva/ocn.sys.mwPerKm
        if (sqrt(ocn.owpps[1].zone.area) < dif)
            dif=sqrt(ocn.owpps[1].zone.area)
        else
        end
        ocn.owpps[1].zone.pos_height=dif/2-ocn.sys.prec/2
        ocn.owpps[1].zone.neg_height=dif/2-ocn.sys.prec/2
        ocn.owpps[1].zone.pos_width=ocn.owpps[1].zone.area/((ocn.owpps[1].zone.pos_height+ocn.owpps[1].zone.neg_height)*2)
        ocn.owpps[1].zone.neg_width=ocn.owpps[1].zone.pos_width

        dif=abs(ocn.owpps[length(ocn.owpps)].node.xy.y-ocn.owpps[length(ocn.owpps)-1].node.xy.y)
        ocn.owpps[length(ocn.owpps)].zone.area=ocn.owpps[length(ocn.owpps)].mva/ocn.sys.mwPerKm
        if (sqrt(ocn.owpps[length(ocn.owpps)].zone.area) < dif)
            dif=sqrt(ocn.owpps[length(ocn.owpps)].zone.area)
        else
        end
        ocn.owpps[length(ocn.owpps)].zone.pos_height=dif/2-ocn.sys.prec/2
        ocn.owpps[length(ocn.owpps)].zone.neg_height=dif/2-ocn.sys.prec/2
        ocn.owpps[length(ocn.owpps)].zone.pos_width=ocn.owpps[length(ocn.owpps)].zone.area/((ocn.owpps[length(ocn.owpps)].zone.pos_height+ocn.owpps[length(ocn.owpps)].zone.neg_height)*2)
        ocn.owpps[length(ocn.owpps)].zone.neg_width=ocn.owpps[length(ocn.owpps)].zone.pos_width

        for indx=2:length(ocn.owpps)-1

            dif=abs(ocn.owpps[indx].node.xy.y-ocn.owpps[indx+1].node.xy.y)
            ocn.owpps[indx].zone.area=ocn.owpps[indx].mva/ocn.sys.mwPerKm
            if (sqrt(ocn.owpps[indx].zone.area) < dif)
                dif=sqrt(ocn.owpps[indx].zone.area)
            else
            end
            #area, height and width
            ocn.owpps[indx].zone.neg_height=ocn.owpps[indx-1].zone.pos_height
            ocn.owpps[indx].zone.pos_height=dif/2-ocn.sys.prec/2
            ocn.owpps[indx].zone.pos_width=ocn.owpps[indx].zone.area/((ocn.owpps[indx].zone.pos_height+ocn.owpps[indx].zone.neg_height)*2)
            ocn.owpps[indx].zone.neg_width=ocn.owpps[indx].zone.pos_width
        end
    else
        #set areas, heights and widths
        dif=abs(ocn.owpps[1].node.xy.x-ocn.owpps[2].node.xy.x)
        ocn.owpps[1].zone.area=ocn.owpps[1].mva/ocn.sys.mwPerKm
        if (sqrt(ocn.owpps[1].zone.area) < dif)
            dif=sqrt(ocn.owpps[1].zone.area)
        else
        end
        ocn.owpps[1].zone.pos_width=dif/2-ocn.sys.prec/2
        ocn.owpps[1].zone.neg_width=dif/2-ocn.sys.prec/2
        ocn.owpps[1].zone.pos_height=ocn.owpps[1].zone.area/((ocn.owpps[1].zone.pos_width+ocn.owpps[1].zone.neg_width)*2)
        ocn.owpps[1].zone.neg_height=ocn.owpps[1].zone.pos_height

        dif=abs(ocn.owpps[length(ocn.owpps)].node.xy.x-ocn.owpps[length(ocn.owpps)-1].node.xy.x)
        ocn.owpps[length(ocn.owpps)].zone.area=ocn.owpps[length(ocn.owpps)].mva/ocn.sys.mwPerKm
        if (sqrt(ocn.owpps[length(ocn.owpps)].zone.area) < dif)
            dif=sqrt(ocn.owpps[length(ocn.owpps)].zone.area)
        else
        end
        ocn.owpps[length(ocn.owpps)].zone.pos_width=dif/2-ocn.sys.prec/2
        ocn.owpps[length(ocn.owpps)].zone.neg_width=dif/2-ocn.sys.prec/2
        ocn.owpps[length(ocn.owpps)].zone.pos_height=ocn.owpps[length(ocn.owpps)].zone.area/((ocn.owpps[length(ocn.owpps)].zone.pos_width+ocn.owpps[length(ocn.owpps)].zone.neg_width)*2)
        ocn.owpps[length(ocn.owpps)].zone.neg_height=ocn.owpps[length(ocn.owpps)].zone.pos_height

        for indx=2:length(ocn.owpps)-1

            dif=abs(ocn.owpps[indx].node.xy.x-ocn.owpps[indx+1].node.xy.x)
            ocn.owpps[indx].zone.area=ocn.owpps[indx].mva/ocn.sys.mwPerKm
            if (sqrt(ocn.owpps[indx].zone.area) < dif)
                dif=sqrt(ocn.owpps[indx].zone.area)
            else
            end
            #area, height and width
            ocn.owpps[indx].zone.neg_width=ocn.owpps[indx-1].zone.pos_width
            ocn.owpps[indx].zone.pos_width=dif/2-ocn.sys.prec/2
            ocn.owpps[indx].zone.pos_height=ocn.owpps[indx].zone.area/((ocn.owpps[indx].zone.pos_width+ocn.owpps[indx].zone.neg_width)*2)
            ocn.owpps[indx].zone.neg_height=ocn.owpps[indx].zone.pos_height
        end
    end
end
#=
ocn=ocean
ngs=ocn.nogos[1]
=#
function lof_bndafy(ocn)
    ocn.nbnd,ocn.ebnd,ocn.sbnd,ocn.wbnd=lof_bndNESW(ocn.bndryPnts)
    for ngs in ocn.nogos
        ngs.nbnd,ngs.ebnd,ngs.sbnd,ngs.wbnd,ngs.bndryPnts=lof_bndNESW_nogo(ngs.bndryPnts)
    end

    ngs_dummy=owppNGnodes(ocn, false)
    for (indx,owpp) in enumerate(ocn.owpps)
        owpp.zone.nbnd,owpp.zone.ebnd,owpp.zone.sbnd,owpp.zone.wbnd,owp_pnts=lof_bndNESW_nogo(ngs_dummy[indx].nodes)
    end

    ngs_dummy=owppNGnodes(ocn, true)
    for (indx,owpp) in enumerate(ocn.owpps)
        owpp.mv_zone.nbnd,owpp.mv_zone.ebnd,owpp.mv_zone.sbnd,owpp.mv_zone.wbnd,owp_pnts=lof_bndNESW_nogo(ngs_dummy[indx].nodes)
    end
end
#=
nds=ocn.bndryPnts
=#
#finds line equations and x/y limits of each for all shapes
function lof_bndNESW(nds)

    vert=Array{line,1}()
    hori=Array{line,1}()
    set_x=Array{Float32,1}()
    set_y=Array{Float32,1}()

    for indx=1:length(nds)
        y0=nds[indx].xy.y
        x0=nds[indx].xy.x
        if (indx != length(nds))
            y1=nds[indx+1].xy.y
            x1=nds[indx+1].xy.x
        else
            y1=nds[1].xy.y
            x1=nds[1].xy.x
        end

        dummy_line=line()
        dummy_line.ymn=min(y0,y1)
        dummy_line.ymx=max(y0,y1)
        dummy_line.xmn=min(x0,x1)
        dummy_line.xmx=max(x0,x1)
        push!(set_x,x0)
        push!(set_x,x1)
        push!(set_y,y0)
        push!(set_y,y1)

        if (abs(y1-y0) > abs(x1-x0))
            #vertical
            #store all vertical xs and ys
            dummy_line.b_findx,dummy_line.m_findx=reverse([[y0,y1] ones(2)]\[x0,x1])
            if (x0==x1)
                dummy_line.b_findy,dummy_line.m_findy=reverse([[x0,x1+(10^-5)] ones(2)]\[y0,y1])
            else
                dummy_line.b_findy,dummy_line.m_findy=reverse([[x0,x1] ones(2)]\[y0,y1])
            end
            push!(vert,dummy_line)
        else
            #horizontal
            println(string(x0)*", "*string(y0))
            println(string(x1)*", "*string(y1))
            #store all horizontal xs and ys
            dummy_line.b_findy,dummy_line.m_findy=reverse([[x0,x1] ones(2)]\[y0,y1])
            if (y0==y1)
                dummy_line.b_findx,dummy_line.m_findx=reverse([[y0,y1+(10^-5)] ones(2)]\[x0,x1])
            else
                dummy_line.b_findx,dummy_line.m_findx=reverse([[y0,y1] ones(2)]\[x0,x1])
            end
            push!(hori,dummy_line)
        end
    end

    #find average x and y of shape
    x_mean=(sum(set_x)/length(set_x))
    y_mean=(sum(set_y)/length(set_y))
    x_max=findmax(set_x)[1]
    y_max=findmax(set_y)[1]
    x_min=findmin(set_x)[1]
    y_min=findmin(set_y)[1]

    #Create arrays to return
    nb=Array{line,1}()
    eb=Array{line,1}()
    sb=Array{line,1}()
    wb=Array{line,1}()

    #sort north and south
    for ln0 in hori
        for lns in [vert,hori]
            for ln1 in lns
                if (ln0 != ln1)
                    for x=ln0.xmn:0.1:ln0.xmx
                        if (x >= ln1.xmn && x <= ln1.xmx)
                            if (x*ln0.m_findy+ln0.b_findy)<(x*ln1.m_findy+ln1.b_findy)
                                push!(sb,ln0)
                                @goto horiz_line_stored
                            end
                        end
                    end
                end
            end
        end
        push!(nb,ln0)
        @label horiz_line_stored
    end

    #sort left and right
    for ln0 in vert
        for lns in [vert,hori]
            for ln1 in lns
                if (ln0 != ln1)
                    for y=ln0.ymn:0.1:ln0.ymx
                        if (y >= ln1.ymn && y <= ln1.ymx)
                            if (y*ln0.m_findx+ln0.b_findx)<(y*ln1.m_findx+ln1.b_findx)
                                push!(wb,ln0)
                                @goto vert_line_stored
                            end
                        end
                    end
                end
            end
        end
        push!(eb,ln0)
        @label vert_line_stored
    end

    return nb,eb,sb,wb
end
#=
lngth_owpps=Array{Tuple{Float64,Int64},1}()
ordrd_owpps=Array{bus,1}()
for (indx,owp) in enumerate(ocn.owpps)
    owp.num=deepcopy(indx)
    push!(lngth_owpps,(deepcopy(lof_pnt2pnt_dist(owp.node.xy,pcc.node.xy)),deepcopy(owp.num)))
end
sort!(lngth_owpps, by = x -> x[1])

ngs=ocean.nogos[1]
nds=ngs.bndryPnts
=#


function lof_bndNESW_nogo(nds)
    #Create arrays to return
    nb=Array{line,1}()
    eb=Array{line,1}()
    sb=Array{line,1}()
    wb=Array{line,1}()
    nds=lof_sortNogo(nds)
    x0=nds[1].xy.x
    x1=nds[2].xy.x
    x2=nds[3].xy.x
    x3=nds[4].xy.x
    y0=nds[1].xy.y
    y1=nds[2].xy.y
    y2=nds[3].xy.y
    y3=nds[4].xy.y
    ang12=(atan((y1-y0)/(x1-x0))*(180/pi))#90>=ang12>45(east)||45>=ang12>=-45(south)||-45>ang12>=-90(west)
    ang23=(atan((y2-y1)/(x2-x1))*(180/pi))#0>ang23>-45(south)||-45>=ang23>=-90(west)||90>=ang23=>45(west)||45>ang23=>0(north)
    ang34=(atan((y3-y2)/(x3-x2))*(180/pi))#90=>ang34>45(east)||45=>ang34(north)=>-45||-45>ang34=>-90(west)
    ang41=(atan((y0-y3)/(x0-x3))*(180/pi))#0=<ang41<45(south)||45=<ang41=<90(east)||0>ang41=>-45(north)||-45>ang41>=-90(east)
#90>=ang12>45(west)||45>=ang12>=-45(south)||-45>ang12>=-90(east)
    if (-91.0<=ang12 && ang12<-45.0)
        push!(eb,lof_makeDline(x0,x1,y0,y1,true))
    elseif (-45.0<=ang12 && ang12<=45.0)
        push!(nb,lof_makeDline(x0,x1,y0,y1,false))
    elseif (45.0<ang12 && ang12<=91.0)
        push!(wb,lof_makeDline(x0,x1,y0,y1,true))
    else
        println("Error: ang12 not in -90 to 90 degree range.")
    end
#0>ang23>-45(south)||-45>=ang23>=-90(west)||90>=ang23=>45(west)||45>ang23=>0(north)
    if (-91.0<=ang23 && ang23<-45.0)
        push!(eb,lof_makeDline(x1,x2,y1,y2,true))
    elseif (-45.0<=ang23 && ang23<0.0)
        push!(nb,lof_makeDline(x1,x2,y1,y2,false))
    elseif (0.0<=ang23 && ang23<45.0)
        push!(sb,lof_makeDline(x1,x2,y1,y2,false))
    elseif (45.0<=ang23 && ang23<=91.0)
        push!(eb,lof_makeDline(x1,x2,y1,y2,true))
    else
        println("Error: ang23 not in -90 to 90 degree range.")
        println(ang23)
    end
#90=>ang34>45(east)||45=>ang34(north)=>-45||-45>ang34=>-90(west)
    if (-91.0<=ang34 && ang34<-45.0)
        push!(wb,lof_makeDline(x2,x3,y2,y3,true))
    elseif (-45.0<=ang34 && ang34<=45.0)
        push!(sb,lof_makeDline(x2,x3,y2,y3,false))
    elseif (45.0<ang34 && ang34<=91.0)
        push!(eb,lof_makeDline(x2,x3,y2,y3,true))
    else
        println("Error: ang34 not in -90 to 90 degree range.")
    end
#0=<ang41<45(south)||45=<ang41=<90(east)||0>ang41=>-45(north)||-45>ang41>=-90(east)
    if (-91.0<=ang41 && ang41<-45.0)
        push!(wb,lof_makeDline(x3,x0,y3,y0,true))
    elseif (-45.0<=ang41 && ang41<0.0)
        push!(sb,lof_makeDline(x3,x0,y3,y0,false))
    elseif (0.0<=ang41 && ang41<45.0)
        push!(nb,lof_makeDline(x3,x0,y3,y0,false))
    elseif (45.0<=ang41 && ang41<=91.0)
        push!(wb,lof_makeDline(x3,x0,y3,y0,true))
    else
        println("Error: ang41 not in -90 to 90 degree range.")
        println(ang41)
    end
    return nb,eb,sb,wb,nds
end

function lof_makeDline(x0,x1,y0,y1,vert)
    dummy_line=line()
    dummy_line.ymn=min(y0,y1)
    dummy_line.ymx=max(y0,y1)
    dummy_line.xmn=min(x0,x1)
    dummy_line.xmx=max(x0,x1)
    if (vert==true)
        dummy_line.b_findx,dummy_line.m_findx=reverse([[y0,y1] ones(2)]\[x0,x1])
        if (x0==x1)
            dummy_line.b_findy,dummy_line.m_findy=reverse([[x0,x1+(10^-5)] ones(2)]\[y0,y1])
        else
            dummy_line.b_findy,dummy_line.m_findy=reverse([[x0,x1] ones(2)]\[y0,y1])
        end
    else
        dummy_line.b_findy,dummy_line.m_findy=reverse([[x0,x1] ones(2)]\[y0,y1])
        if (y0==y1)
            dummy_line.b_findx,dummy_line.m_findx=reverse([[y0,y1+(10^-5)] ones(2)]\[x0,x1])
        else
            dummy_line.b_findx,dummy_line.m_findx=reverse([[y0,y1] ones(2)]\[x0,x1])
        end
    end
    return dummy_line
end

function lof_sortNogo(pnts)
    _nodes=deepcopy(pnts)
    sorted_nodes=node[]
    _nodeT=Array{Tuple{node,Int64},1}()
    n_xy=Array{Tuple{xy,Int64},1}()
    s_xy=Array{Tuple{xy,Int64},1}()
    _xy=Array{Tuple{xy,Int64},1}()
    sorted_xy=Array{Tuple{xy,Int64},1}()
    _x=Array{Tuple{Float64,Int64},1}()
    _y=Array{Tuple{Float64,Int64},1}()
    for (i,ns) in enumerate(_nodes)
        push!(_nodeT,(ns,i))
        push!(_xy,(ns.xy,i))
        push!(_y,(ns.xy.y,i))
    end
    sort!(_y, by = x -> x[1])
    for y in _y
        for xys in _xy
            if (xys[2]==y[2])
                push!(sorted_xy,xys)
            end
        end
    end
    s_xy=sorted_xy[1:2]
    n_xy=sorted_xy[3:4]
    if (n_xy[1][1].x>n_xy[2][1].x)
        n_xy=n_xy[end:-1:1]
    end
    if (s_xy[1][1].x<s_xy[2][1].x)
        s_xy=s_xy[end:-1:1]
    end
    for nxy in n_xy
        for nd in _nodeT
            if (nd[2]==nxy[2])
                push!(sorted_nodes,nd[1])
            end
        end
    end
    for sxy in s_xy
        for nd in _nodeT
            if (nd[2]==sxy[2])
                push!(sorted_nodes,nd[1])
            end
        end
    end
    return sorted_nodes
end

function lof_nodify(ocn)
    wbnd,ebnd=lof_bndPerimeter(ocn.bndryPnts)
    Allnodes=lof_addNodes(ocn,wbnd,ebnd)
    Allnodes=lof_deleteNgNodes(ocn.nogos,Allnodes)
    ngs_dummy=owppNGnodes(ocn)
    ocn.discretedom.nodes=lof_deleteNgNodes(ngs_dummy,Allnodes)


    #lof_busPlaceOnNodes(ocn.owpps,ocn.discretedom.nodes)
    #lof_busPlaceOnNodes(ocn.pccs,ocn.discretedom.nodes)
    #lof_nogoAreaNodes(ocn)
    lof_owppAreaNodes(ocn)
    lof_mvAreaNodes(ocn)
end

function owppNGnodes(ocn, mv)
    ngs_dummy=Array{nogo,1}()
    for owpp in ocn.owpps
        if mv==true
            zone=owpp.mv_zone
        else
            zone=owpp.zone
        end
        nogo_dummy=nogo()
        node_dummy=node()

        node_dummy.xy.y=owpp.node.xy.y-zone.neg_height
        node_dummy.xy.x=owpp.node.xy.x-zone.neg_width
        push!(nogo_dummy.nodes,deepcopy(node_dummy))

        node_dummy.xy.y=owpp.node.xy.y-zone.neg_height
        node_dummy.xy.x=owpp.node.xy.x+zone.pos_width
        push!(nogo_dummy.nodes,deepcopy(node_dummy))

        node_dummy.xy.y=owpp.node.xy.y+zone.pos_height
        node_dummy.xy.x=owpp.node.xy.x+zone.pos_width
        push!(nogo_dummy.nodes,deepcopy(node_dummy))

        node_dummy.xy.y=owpp.node.xy.y+zone.pos_height
        node_dummy.xy.x=owpp.node.xy.x-zone.neg_width
        push!(nogo_dummy.nodes,deepcopy(node_dummy))

        push!(ngs_dummy,deepcopy(nogo_dummy))
    end

    return ngs_dummy
end

#=
function lof_nogoPerimeter(ngs)
    wbndNG=Array{line,1}()
    ebndNG=Array{line,1}()
    for i=1:length(ngs)
        wbnd_ng,ebnd_ng=lof_bndPerimeter(ngs[i].nodes)
        for ln in wbnd_ng
            push!(wbndNG,deepcopy(ln))
        end
        for ln in ebnd_ng
            push!(ebndNG,deepcopy(ln))
        end
    end
    return wbndNG, ebndNG
end
=#

function lof_nogoAreaNodes(ocn)
    for node in ocn.discretedom.nodes
        for owpp in ocn.owpps
            km2owpp=lof_pnt2pnt_dist(node.xy,owpp.node.xy)
            if km2owpp<=owpp.zone.radius
                push!(owpp.zone.pnts,node)
                if km2owpp+ocn.sys.prec>owpp.zone.radius
                    push!(owpp.zone.pnts,node)
                end
            end
        end
    end
end

#=
function lof_busPlaceOnNodes(buses,nodes)
    dummy_node=node()
    for bus in buses
        closestNode=Inf
        for node in nodes
            km2owpp=lof_pnt2pnt_dist(node.xy,bus.node.xy)
            if km2owpp<closestNode
                closestNode=deepcopy(km2owpp)
                dummy_node=node
            end
        end
        bus.node=deepcopy(dummy_node)
    end
end
=#

function lof_owppAreaNodes(ocn)
    for owpp in ocn.owpps

        const_sy=owpp.node.xy.y-owpp.zone.neg_height
        const_ny=owpp.node.xy.y+owpp.zone.pos_height
        mn_y=const_sy
        mx_y=const_ny

        const_wx=owpp.node.xy.x-owpp.zone.neg_width
        const_ex=owpp.node.xy.x+owpp.zone.pos_width
        mn_x=const_wx
        mx_x=const_ex

        for indy=mn_y:ocn.sys.prec/2:mx_y
            dummy_wnode=Int8
            dummy_enode=Int8
            wbsf_km=Inf
            ebsf_km=Inf
            for (val, node) in enumerate(ocn.discretedom.nodes)
                km2w=lof_pnt2pnt_dist(node.xy,xy(const_wx,indy))
                km2e=lof_pnt2pnt_dist(node.xy,xy(const_ex,indy))
                if (km2e<ebsf_km)
                    ebsf_km=deepcopy(km2e)
                    dummy_enode=deepcopy(val)
                end
                if (km2w<wbsf_km)
                    wbsf_km=deepcopy(km2w)
                    dummy_wnode=deepcopy(val)
                end
            end
            push!(owpp.zone.pnts,ocn.discretedom.nodes[dummy_enode])
            push!(owpp.zone.pnts,ocn.discretedom.nodes[dummy_wnode])
        end
        for indx=mn_x:ocn.sys.prec/2:mx_x
            dummy_nnode=Int8
            dummy_snode=Int8
            nbsf_km=Inf
            sbsf_km=Inf
            for (val,node) in enumerate(ocn.discretedom.nodes)
                km2n=lof_pnt2pnt_dist(node.xy,xy(indx,const_ny))
                km2s=lof_pnt2pnt_dist(node.xy,xy(indx,const_sy))
                if (km2n<nbsf_km)
                    nbsf_km=deepcopy(km2n)
                    dummy_nnode=deepcopy(val)
                end
                if (km2s<sbsf_km)
                    sbsf_km=deepcopy(km2s)
                    dummy_snode=deepcopy(val)
                end
            end
            push!(owpp.zone.pnts,ocn.discretedom.nodes[dummy_snode])
            push!(owpp.zone.pnts,ocn.discretedom.nodes[dummy_nnode])
        end
        unique!(owpp.zone.pnts)
    end
end

function lof_mvAreaNodes(ocn)
    for node in ocn.discretedom.nodes
        for owpp in ocn.owpps
            s=owpp.node.xy.y-owpp.mv_zone.neg_height
            n=owpp.node.xy.y+owpp.mv_zone.pos_height
            w=owpp.node.xy.x-owpp.mv_zone.neg_width
            e=owpp.node.xy.x+owpp.mv_zone.pos_width
            if ((node.xy.y<n) && (node.xy.y>s) && (node.xy.x < e) && (node.xy.x > w))
                push!(owpp.mv_zone.pnts,node)
            end
        end
    end
end
#=
nogoReg=ngs_dummy
ngs=nogoReg[2]
=#
function lof_deleteNgNodes(nogoReg,Allnodes)
    for ngs in nogoReg
        wbnd,ebnd=lof_bndPerimeter(ngs.nodes)
        ymx=wbnd[length(wbnd)].ymax
        ymn=wbnd[1].ymn
        windex=1
        eindex=1
        lnth=length(Allnodes)
        i=1
        while i<lnth

            if ((Allnodes[i].xy.y < ymx) && (Allnodes[i].xy.y > ymn))
                while (Allnodes[i].xy.y >= wbnd[windex].ymax)
                    windex=windex+1
                end
                while (Allnodes[i].xy.y >= ebnd[eindex].ymax)
                    eindex=eindex+1
                end
                xmx=Allnodes[i].xy.y*ebnd[eindex].m+ebnd[eindex].b
                xmn=Allnodes[i].xy.y*wbnd[windex].m+wbnd[windex].b
                if ((Allnodes[i].xy.x > xmn) && (Allnodes[i].xy.x < xmx))
                    deleteat!(Allnodes,i)
                    i=i-1
                    lnth=lnth-1
                else
                end
            else
            end
            i=i+1
        end
    end
    return deepcopy(Allnodes)
end

function lof_addNodes(ocn,wbnd,ebnd)
    dummy_nodes=Array{node,1}()
    ymx=wbnd[length(wbnd)].ymax
    ymn=wbnd[1].ymn
    prec=ocn.sys.prec
    windex=1
    eindex=1
    for y=ymn:prec:ymx
        while y>wbnd[windex].ymax
            windex=windex+1
        end
        while y>ebnd[eindex].ymax
            eindex=eindex+1
        end
        xmx=y*ebnd[eindex].m+ebnd[eindex].b
        xmn=y*wbnd[windex].m+wbnd[windex].b
        x=xmn
        while (x>=xmn && x<=xmx)
            dummy_node=node()
            dummy_node.xy.x=x
            dummy_node.xy.y=y
            push!(dummy_nodes,deepcopy(dummy_node))
            x=x+prec
        end
    end
    return dummy_nodes
end

#=
bnd=ngs.nodes
bnd=ocn.bndryPnts
=#
function lof_bndPerimeter(bnd)
    #finds mid line dividing east and west regions

    ys=Array{Float32,1}()
    for xys in bnd
        push!(ys,xys.xy.y)
    end
    mxY=findmax(ys)
    mnY=findmin(ys)
    midLn=reverse([[bnd[mxY[2]].xy.y,bnd[mnY[2]].xy.y] ones(2)]\[bnd[mxY[2]].xy.x,bnd[mnY[2]].xy.x])

    #sorts to east and west bounding points
    wpntsT=Array{xy,1}()
    epntsT=Array{xy,1}()
    wys=Array{Float32,1}()
    eys=Array{Float32,1}()
    wxs=Array{Float32,1}()
    exs=Array{Float32,1}()
    push!(wpntsT,bnd[mnY[2]].xy)
    push!(epntsT,bnd[mnY[2]].xy)
    push!(wxs, bnd[mnY[2]].xy.x)
    push!(exs, bnd[mnY[2]].xy.x)
    for xys in bnd
        if xys.xy.x < (xys.xy.y*midLn[2]+midLn[1])
            push!(wpntsT, xys.xy)
            push!(wxs, xys.xy.x)
        elseif xys.xy.x > (xys.xy.y*midLn[2]+midLn[1])
            push!(epntsT, xys.xy)
            push!(exs, xys.xy.x)
        else
        end
    end
    push!(wpntsT,bnd[mxY[2]].xy)
    push!(epntsT,bnd[mxY[2]].xy)
    push!(wxs, bnd[mxY[2]].xy.x)
    push!(exs, bnd[mxY[2]].xy.x)

    #orders points xmin to xmax
    xwpnts=Array{xy,1}()
    xepnts=Array{xy,1}()
    for i =1:length(wxs)
        index=findmin(wxs)[2]
        wxs[index]=Inf
        push!(xwpnts,wpntsT[index])
        push!(wys,wpntsT[index].y)
    end
    for i =1:length(exs)
        index=findmin(exs)[2]
        exs[index]=Inf
        push!(xepnts,epntsT[index])
        push!(eys,epntsT[index].y)
    end

    #orders points ymin to ymax
    wpnts=Array{xy,1}()
    epnts=Array{xy,1}()
    for i =1:length(wys)
        index=findmin(wys)[2]
        wys[index]=Inf
        push!(wpnts,xwpnts[index])
    end
    for i =1:length(eys)
        index=findmin(eys)[2]
        eys[index]=Inf
        push!(epnts,xepnts[index])
    end

    #Keep only unique entries
    unique!(epnts)
    unique!(wpnts)

    #constructs line boundaries for limits
    wbnd=Array{line,1}()
    ebnd=Array{line,1}()
    for i = 1:length(wpnts)-1
        dummy_line=line()
        if (wpnts[i+1].y != wpnts[i].y)
            alpha_beta=reverse([[wpnts[i+1].y,wpnts[i].y] ones(2)]\[wpnts[i+1].x,wpnts[i].x])
            dummy_line.b=alpha_beta[1]
            dummy_line.m=alpha_beta[2]
        else
            dummy_line.b=wpnts[i].y
            dummy_line.m=0
        end
        dummy_line.ymax=wpnts[i+1].y
        dummy_line.ymn=wpnts[i].y
        push!(wbnd,deepcopy(dummy_line))
    end
    for i = 1:length(epnts)-1
        dummy_line=line()
        if (epnts[i+1].y != epnts[i].y)
            alpha_beta=reverse([[epnts[i+1].y,epnts[i].y] ones(2)]\[epnts[i+1].x,epnts[i].x])
            dummy_line.b=alpha_beta[1]
            dummy_line.m=alpha_beta[2]
        else
            dummy_line.b=epnts[i].y
            dummy_line.m=0
        end
        dummy_line.ymax=epnts[i+1].y
        dummy_line.ymn=epnts[i].y
        push!(ebnd,deepcopy(dummy_line))
    end
    return wbnd,ebnd
end

function lof_getPccData(id_count)
    df = DataFrame(XLSX.readtable("layout//data.xlsx", "pcc_data")...)
    pccs=Array{bus,1}()
    for index=1:length(df[:, 1])
        dummy_bus=bus()
        dummy_bus.node.gps.lng=df.longitude[index]
        dummy_bus.node.gps.lat=df.latitude[index]
        dummy_bus.kV=df.kv[index]
        dummy_bus.id=id_count
        push!(pccs,deepcopy(dummy_bus))
        pccs[length(pccs)].num=length(pccs)
        id_count=id_count+1
    end
    return pccs,id_count
end

function lof_getOwppData(id_count)
    df = DataFrame(XLSX.readtable("layout//data.xlsx", "owpp_data")...)
    owpps=Array{bus,1}()
    for index=1:length(df[:, 1])
        dummy_bus=bus()
        dummy_bus.node.gps.lng=df.longitude[index]
        dummy_bus.node.gps.lat=df.latitude[index]
        dummy_bus.name=df.name[index]
        dummy_bus.id=id_count
        dummy_bus.mva=df.mva[index]
        push!(owpps,deepcopy(dummy_bus))
        owpps[length(owpps)].num=length(owpps)
        id_count=id_count+1
    end
    lof_getOwppWindData(owpps)
    return owpps,id_count
end

function lof_getOwppWindData(owpps)
    df = DataFrame(XLSX.readtable("layout//data.xlsx", "wind_data")...)
    for owpp in owpps
        wnd=wndF_wndPrf([getproperty(df,Symbol(owpp.name))])
        owpp.wnd=deepcopy(wnd)
    end
end

function lof_getBndData()
    df = DataFrame(XLSX.readtable("layout//data.xlsx", "domain_data")...)
    boundary=Array{node,1}()
    for index=1:length(df[:, 1])
        dummy_node=node()
        dummy_node.gps.lng=df.longitude[index]
        dummy_node.gps.lat=df.latitude[index]
        push!(boundary,deepcopy(dummy_node))
    end
    return boundary
end

function lof_getNoGoData(index)
    df = DataFrame(XLSX.readtable("layout//data.xlsx", "nogo_data"*string(index))...)
    nogo=Array{node,1}()
    for index=1:length(df[:, 1])
        dummy_node=node()
        dummy_node.gps.lng=df.longitude[index]
        dummy_node.gps.lat=df.latitude[index]
        push!(nogo,deepcopy(dummy_node))
    end
    return nogo
end

function lof_getSysData()
    sht = XLSX.readxlsx("layout//data.xlsx")["sys_data"]
    sys=system()
    sys.nogoNum=sht["B1"]
    sys.prec=sht["B2"]
    sys.mwPerKm=sht["B3"]
    return sys
end
############################ GPS to cartesian transform ########################
################################################################################
#sets the gps coords that used as reference coords
function lof_bseCrd(ocean)
    base=gps()
    #for india (type layouts)
    if ocean.owpps[length(ocean.owpps)].node.gps.lat < ocean.pccs[length(ocean.pccs)].node.gps.lat
        base.lat=ocean.owpps[length(ocean.owpps)].node.gps.lat#base lat
        base.lng=ocean.owpps[length(ocean.owpps)].node.gps.lng#base long
    #for belgium (type layouts)
elseif ocean.pccs[length(ocean.pccs)].node.gps.lat < ocean.owpps[length(ocean.owpps)].node.gps.lat
        base.lat=ocean.pccs[length(ocean.pccs)].node.gps.lat#base lat
        base.lng=ocean.pccs[length(ocean.pccs)].node.gps.lng#base long
    else
        error("No proper base coordinates system established!")
    end
    return base
end

#calculates lengths based on latitude
#as lattitude changes number of km should be updated
function lof_gps2cartesian4nodes(location,base)
    lnthLT=111#number of km in 1 degree of longitude at equator
    for value in location
        value.xy.x=lof_deg2lgth(value.gps.lng-base.lng,lof_lg1deg(value.gps.lat,lnthLT))
        value.xy.y=lof_deg2lgth(value.gps.lat-base.lat,lnthLT)
    end
end

#as lattitude changes number of km should be updated
function lof_gps2cartesian4bus(location,base)
    lnthLT=111#number of km in 1 degree of longitude at equator
    for value in location
        value.node.xy.x=lof_deg2lgth(value.node.gps.lng-base.lng,lof_lg1deg(value.node.gps.lat,lnthLT))
        value.node.xy.y=lof_deg2lgth(value.node.gps.lat-base.lat,lnthLT)
    end
end

#rotates and slides cartesian axis
function lof_transformAxis(ocn)
    offset=lof_rotateAxis(ocn)
    lof_slideAxis(ocn,offset)
end

#finds angle to rotate and applies to owpps and pccs
#rotates axis to align n-s with y
function lof_rotateAxis(ocn)
    theta=atan((ocn.pccs[length(ocn.pccs)].node.xy.x-ocn.owpps[length(ocn.owpps)].node.xy.x)/(ocn.owpps[length(ocn.owpps)].node.xy.y-ocn.pccs[length(ocn.pccs)].node.xy.y))
    ocn.theta=theta
    offset=0.0
    offset=lof_rotateGroup4bus(ocn.owpps,theta,offset)
    offset=lof_rotateGroup4bus(ocn.pccs,theta,offset)
    ocn.offset=deepcopy(offset)
    lof_rotateGroup4node(ocn.bndryPnts,theta)
    for i=1:length(ocn.nogos)
        lof_rotateGroup4node(ocn.nogos[i].bndryPnts,theta)
    end
    return ocn.offset
end

#loops through to apply rotations for a specified group
function lof_rotateGroup4bus(locations,theta,os)
    for value in locations
        xy=lof_rotatePnt(value.node.xy.x,value.node.xy.y,theta)
        value.node.xy.x=xy[1]
        value.node.xy.y=xy[2]
        if value.node.xy.x<os
            os=value.node.xy.x
        end
    end
    return os
end

#loops through to apply rotations for a specified group
function lof_rotateGroup4node(locations,theta)
    for value in locations
        xy=lof_rotatePnt(value.xy.x,value.xy.y,theta)
        value.xy.x=xy[1]
        value.xy.y=xy[2]
    end
end

#applies rotational matrix individual coordinates
function lof_rotatePnt(x,y,theta)
    co_od=[x y]
    rotated=co_od*[cos(theta) -1*sin(theta);sin(theta) cos(theta)]
    return rotated
end

#translates the entire region by specified offset
#sets unique IDs for owpps and pccs
function lof_slideAxis(ocn,os)
    for value in ocn.owpps
        value.node.xy.x=value.node.xy.x-os
    end
    for value in ocn.pccs
        value.node.xy.x=value.node.xy.x-os
    end
    for value in ocn.bndryPnts
        value.xy.x=value.xy.x-os
    end
    for i=1:length(ocn.nogos)
        for value in ocn.nogos[i].bndryPnts
            value.xy.x=value.xy.x-os
        end
    end
end

#changes angle to an arc length
function lof_deg2lgth(d,dPl)
    return d*dPl
end

#calculates length of 1 deg of longitude at given lattitude
function lof_lg1deg(lat,lngth)
    return cos(lof_d2r(lat))*lngth
end

#Change radians to degrees
function lof_r2d(rad)
    return rad*180/pi
end

#Change degrees to radians
function lof_d2r(deg)
    return deg*pi/180
end


#returns the hypotenuse distance between 2 cartesian points
#minimum distance for a path is 1km
function lof_pnt2pnt_dist(pnt1,pnt2)
    hyp=sqrt((pnt2.x-pnt1.x)^2+(pnt2.y-pnt1.y)^2)
    return hyp
end
