
#=
ocn=ocean
owpps=ocn.owpps
pcc=ocn.pccs[2]
paths=opt_mvOSSplacement(ocn,owpps,pcc)
=#
function opt_hvOSSplacement(ocn,pcc)
    mv_circs=ocn.circuits
    owpps_tbl_hv=Array{Array{Int8, 1}, 1}()
    oss_system=circuit()
    oss_systems=Array{circuit,1}()
    for crc in mv_circs
        #if (crc.base_owp.num!=findfirst(x->x==true,crc.binary))
            push!(owpps_tbl_hv,crc.binary)
        #end
    end
    for indx0=1:length(owpps_tbl_hv)
    #for indx0=225:225
        bus_dummies=bus[]
        for (indx1,owp) in enumerate(owpps_tbl_hv[indx0])
            if (owp==1)
                push!(bus_dummies,ocn.owpps[indx1])
            end
        end
        oss_system=opt_hvOssSystem(bus_dummies,pcc,ocn,owpps_tbl_hv[indx0])
        push!(oss_systems,deepcopy(oss_system))
    end
    return oss_systems
end
#length(owpps_tbl[:,1])
#length(owpps_tbl_hv[:,1])
function opt_mvOSSplacement(ocn,owpps,pcc)
    #owpps is an array of owpps ordered closest to farthest from the designated pcc
    owpps_tbl=top_mvTopos(owpps)
    bus_dummies=Array{bus,1}()
    oss_system=circuit()
    oss_systems=Array{circuit,1}()
    for indx0=1:length(owpps_tbl[:,1])
    #for indx0=225:225
        bus_dummies=bus[]
        mv_square=opt_mvConstraints(ocn,owpps_tbl[indx0,:])
        for (indx1,owpp) in enumerate(owpps_tbl[indx0,:])
            if (owpp==1)
                push!(bus_dummies,owpps[indx1])
            end
        end
        if (sum(owpps_tbl[indx0,:])==1)
            oss_system=opt_str8Connect(bus_dummies[1],pcc,ocn,owpps_tbl[indx0,:])
        else
            oss_system=opt_mvOssSystem(mv_square,bus_dummies,pcc,ocn,owpps_tbl[indx0,:])
        end
        push!(oss_systems,deepcopy(oss_system))
    end
    return oss_systems
end
#=
function opt_mvlnth(syss)
    for sys in syss
        for pth in sys.pths
            push!(sys.lengths,pth.G_cost)
        end
    end
end=#
#owp=bus_dummies[1]
#bn=owpps_tbl[indx0,:]
function opt_str8Connect(owp,pcc,ocn,bn)
    circ=circuit()
    circ.binary=bn
    circ.decimal=top_bin2dec(bn)
    circ.pcc=pcc
    push!(circ.owpps,owp)
    push!(circ.pths,as_Astar(owp.node,pcc.node,ocn.discretedom.nodes))
    push!(circ.lengths,circ.pths[1].G_cost)
    #not necessarily MV
    cs=cstF_MvHvCblpcc(circ.lengths[length(circ.lengths)],owp.mva,owp.wnd,ocn.finance,pcc)
    for c in cs
        push!(circ.cbls,c)
    end
    circ.base_owp=owp

    if (length(circ.cbls) > 1)
        push!(circ.xfmrs,cstF_xfo_oss(owp.mva,owp.wnd,ocean.finance))
    end
    if (circ.cbls[1].elec.volt != pcc.kV)
        push!(circ.xfmrs,cstF_xfo_pcc(owp.mva,owp.wnd,ocean.finance))
    end
    opt_ttlMvCirc(circ)
    return circ
end
#buses=bus_dummies
#bn=owpps_tbl_hv[indx0]

function opt_hvOssSystem(buses,pcc,ocn,bn)
    xys=opt_mvhvOss1stLocal(buses,pcc,ocn)
    owp=buses[1]
    domain_oss,oss_node,power_sum,wind_sum=opt_adjustPath(owp,ocn,xys,buses,pcc)
    circ=circuit()
    push!(circ.pths,deepcopy(as_Astar(domain_oss[pcc.node.num],oss_node,domain_oss)))
    push!(circ.lengths,circ.pths[length(circ.pths)].G_cost)
    println("mva: "*string(power_sum))
    cbl_xfo=cstF_HvCblallKvo2p(circ.lengths[length(circ.lengths)],power_sum,wind_sum,ocn.finance,pcc)
    push!(circ.cbls,cbl_xfo[1])
    pSum_hv=0
    pSum_hv132=0
    pSum_hv220=0
    pSum_hv400=0
    pSum_mv=0
    iSum_hv=0
    iSum_hv132=0
    iSum_hv220=0
    iSum_hv400=0
    iSum_mv=0
    wind_sumMv=wind()
    wind_sumMv.pu=zeros(Float64,length(buses[1].wnd.pu))
    wind_sumMv.ce=zeros(Float64,length(buses[1].wnd.ce))
    wind_sumMv.delta=0
    wind_sumMv.lf=0
    wind_sumHv=deepcopy(wind_sumMv)
    wind_sum132=deepcopy(wind_sumMv)
    wind_sum220=deepcopy(wind_sumMv)
    wind_sum400=deepcopy(wind_sumMv)
    for owp in buses
        push!(circ.pths,deepcopy(as_Astar(domain_oss[owp.node.num],oss_node,domain_oss)))
        push!(circ.lengths,circ.pths[length(circ.pths)].G_cost)
        cs,xs=cstF_MvHvCbloss(circ.lengths[length(circ.lengths)],owp.mva,owp.wnd,ocn.finance,pcc,cbl_xfo[1].elec.volt)
        for c in cs
            push!(circ.cbls,c)
        end
        if length(cs) == 1
            iSum_mv=iSum_mv+1
            pSum_mv=pSum_mv+deepcopy(owp.mva)
            wind_sumMv.pu=deepcopy((wind_sumMv.pu.+owp.wnd.pu))
            wind_sumMv.ce=deepcopy((wind_sumMv.ce.+owp.wnd.ce))
            wind_sumMv.delta=deepcopy((wind_sumMv.delta+owp.wnd.delta))
            wind_sumMv.lf=deepcopy((wind_sumMv.lf+owp.wnd.lf))
        elseif cs[2].elec.volt == cbl_xfo[1].elec.volt
            push!(circ.xfmrs,xs[1])
            iSum_hv=iSum_hv+1
            pSum_hv=pSum_hv+10
            wind_sumHv.pu=(wind_sumHv.pu.+deepcopy(owp.wnd.pu))
            wind_sumHv.ce=(wind_sumHv.ce.+deepcopy(owp.wnd.ce))
            wind_sumHv.delta=(wind_sumHv.delta+deepcopy(owp.wnd.delta))
            wind_sumHv.lf=(wind_sumHv.lf+deepcopy(owp.wnd.lf))
        elseif cs[2].elec.volt == 132
            push!(circ.xfmrs,xs[1])
            iSum_hv132=iSum_hv132+1
            pSum_hv132=pSum_hv132+deepcopy(owp.mva)
            wind_sum132.pu=(wind_sum132.pu.+deepcopy(owp.wnd.pu))
            wind_sum132.ce=(wind_sum132.ce.+deepcopy(owp.wnd.ce))
            wind_sum132.delta=(wind_sum132.delta+deepcopy(owp.wnd.delta))
            wind_sum132.lf=(wind_sum132.lf+deepcopy(owp.wnd.lf))
        elseif cs[2].elec.volt == 220
            push!(circ.xfmrs,xs[1])
            iSum_hv220=iSum_hv220+1
            pSum_hv220=pSum_hv220+deepcopy(owp.mva)
            wind_sum220.pu=(wind_sum220.pu.+deepcopy(owp.wnd.pu))
            wind_sum220.ce=(wind_sum220.ce.+deepcopy(owp.wnd.ce))
            wind_sum220.delta=(wind_sum220.delta+deepcopy(owp.wnd.delta))
            wind_sum220.lf=(wind_sum220.lf+deepcopy(owp.wnd.lf))
        elseif cs[2].elec.volt == 400
            push!(circ.xfmrs,xs[1])
            iSum_hv400=iSum_hv400+1
            pSum_hv400=pSum_hv400+deepcopy(owp.mva)
            wind_sum400.pu=(wind_sum400.pu.+deepcopy(owp.wnd.pu))
            wind_sum400.ce=(wind_sum400.ce.+deepcopy(owp.wnd.ce))
            wind_sum400.delta=(wind_sum400.delta+deepcopy(owp.wnd.delta))
            wind_sum400.lf=(wind_sum400.lf+deepcopy(owp.wnd.lf))
        end
    end
    wind_sumMv.pu=wind_sumMv.pu./iSum_mv
    wind_sumMv.ce=wind_sumMv.ce./iSum_mv
    wind_sumMv.delta=wind_sumMv.delta/iSum_mv
    wind_sumMv.lf=wind_sumMv.lf/iSum_mv

    wind_sumHv.pu=wind_sumHv.pu./iSum_hv
    wind_sumHv.ce=wind_sumHv.ce./iSum_hv
    wind_sumHv.delta=wind_sumHv.delta/iSum_hv
    wind_sumHv.lf=wind_sumHv.lf/iSum_hv

    wind_sum132.pu=wind_sum132.pu./iSum_hv132
    wind_sum132.ce=wind_sum132.ce./iSum_hv132
    wind_sum132.delta=wind_sum132.delta/iSum_hv132
    wind_sum132.lf=wind_sum132.lf/iSum_hv132

    wind_sum220.pu=wind_sum220.pu./iSum_hv220
    wind_sum220.ce=wind_sum220.ce./iSum_hv220
    wind_sum220.delta=wind_sum220.delta/iSum_hv220
    wind_sum220.lf=wind_sum220.lf/iSum_hv220

    wind_sum400.pu=wind_sum400.pu./iSum_hv400
    wind_sum400.ce=wind_sum400.ce./iSum_hv400
    wind_sum400.delta=wind_sum400.delta/iSum_hv400
    wind_sum400.lf=wind_sum400.lf/iSum_hv400

    if (pSum_hv132 == 0 && pSum_hv220 == 0 && pSum_hv400 == 0 && pSum_mv == 0)
        pSum_hv=pSum_hv+100
        push!(circ.xfmrs,cstF_xfo_oss(pSum_hv,wind_sumHv,ocean.finance))
    elseif pSum_hv != 0
        push!(circ.xfmrs,cstF_xfo_oss(pSum_hv,wind_sumHv,ocean.finance))
    end

    if pSum_mv != 0
        push!(circ.xfmrs,cstF_xfo_oss(pSum_mv,wind_sumMv,ocean.finance))
    end

    if pSum_hv132 != 0
        push!(circ.xfmrs,cstF_xfo_oss(pSum_hv132,wind_sum132,ocean.finance))
    end

    if pSum_hv220 != 0
        push!(circ.xfmrs,cstF_xfo_oss(pSum_hv220,wind_sum220,ocean.finance))
    end

    if pSum_hv400 != 0
        push!(circ.xfmrs,cstF_xfo_oss(pSum_hv400,wind_sum400,ocean.finance))
    end

    if (pcc.kV != circ.cbls[1].elec.volt)
        push!(circ.xfmrs,cstF_xfo_pcc(power_sum,wind_sum,ocean.finance))
    end


    circ.binary=bn
    circ.decimal=top_bin2dec(bn)
    circ.pcc=pcc
    circ.owpps=buses
    circ.base_owp=owp
    opt_ttlMvCirc(circ)
    return circ
end

#bn=owpps_tbl[indx0,:]

function opt_mvOssSystem(mv_square,buses,pcc,ocn,bn)
    xys=opt_mvhvOss1stLocal(buses,pcc,ocn)
    owp=opt_buses2nodes4MVoptLocal(mv_square,ocn.owpps,buses)
    domain_oss,oss_node,power_sum,wind_sum=opt_adjustPath(owp,ocn,xys,buses,pcc)
    circ=circuit()
    push!(circ.pths,deepcopy(as_Astar(domain_oss[pcc.node.num],oss_node,domain_oss)))
    push!(circ.lengths,circ.pths[length(circ.pths)].G_cost)
    #################### This needs to be a high voltage cable!!!
    println("mva: "*string(power_sum))
    cbl_xfo=cstF_HvCblallKvo2p(circ.lengths[length(circ.lengths)],power_sum,wind_sum,ocn.finance,pcc)
    push!(circ.cbls,cbl_xfo[1])
    push!(circ.xfmrs,cstF_xfo_oss(power_sum,wind_sum,ocean.finance))
    if (circ.cbls[length(circ.cbls)].elec.volt != pcc.kV)
        push!(circ.xfmrs,cbl_xfo[2])
    end
    for owp in buses
        push!(circ.pths,deepcopy(as_Astar(domain_oss[owp.node.num],oss_node,domain_oss)))
        push!(circ.lengths,circ.pths[length(circ.pths)].G_cost)
        push!(circ.cbls,cstF_MvCbl3366(circ.lengths[length(circ.lengths)],owp.mva,owp.wnd,ocn.finance))
    end


    circ.binary=bn
    circ.decimal=top_bin2dec(bn)
    circ.pcc=pcc
    #circ.pcc_length=69.69
    circ.owpps=buses
    #circ.owp_lengths=69.69
    circ.base_owp=owp
    opt_ttlMvCirc(circ)
    return circ
end

function opt_adjustPath(owp,ocn,xys,buses,pcc)
    paths=node[]
    oss_node=opt_OssOptimalLocal(xys,ocn.constrain.ellipses,owp,length(ocn.owpps))
    domain_oss=deepcopy(ocn.discretedom.nodes)
    edges_oss=deepcopy(ocn.discretedom.edges)
    bs,oss_node=opt_ifIn2Out(oss_node,ocn)
    if bs==true
        oss_node.num=deepcopy(length(domain_oss)+1)
        push!(domain_oss,oss_node)
        lof_edgeifyOss(oss_node,ocn,domain_oss,edges_oss)
    end
    push!(paths,deepcopy(as_Astar(domain_oss[pcc.node.num],oss_node,domain_oss)))
    for owpp in buses
        push!(paths,deepcopy(as_Astar(domain_oss[owpp.node.num],oss_node,domain_oss)))
    end

    #pick first node in each owpp path
    xys=xy[]
    pcc_xy=xy()
    #if (paths[1].parent.num != pcc.node.num)
    #    pcc_xy=paths[1].parent.xy
    #else
        pcc_xy=pcc.node.xy
    #end
    push!(xys,pcc_xy)
    for path=2:length(paths)
        pathsParent=paths[path]
        while pathsParent.parent.num != buses[path-1].node.num
            pathsParent=pathsParent.parent
        end
        push!(xys,deepcopy(pathsParent.xy))
    end

    #re calculate the path
    paths=node[]
    oss_node=opt_OssOptimalLocal(xys,ocn.constrain.ellipses,owp,length(ocn.owpps))
    domain_oss=deepcopy(ocn.discretedom.nodes)
    edges_oss=deepcopy(ocn.discretedom.edges)
    bs,oss_node=opt_ifIn2Out(oss_node,ocn)
    if bs==true
        oss_node.num=deepcopy(length(domain_oss)+1)
        push!(domain_oss,oss_node)
        lof_edgeifyOss(oss_node,ocn,domain_oss,edges_oss)
    end
    push!(paths,deepcopy(as_Astar(domain_oss[pcc.node.num],oss_node,domain_oss)))
    power_sum=0
    wind_sum=wind()
    wind_sum.pu=zeros(Float64,length(buses[1].wnd.pu))
    wind_sum.ce=zeros(Float64,length(buses[1].wnd.ce))
    wind_sum.delta=0
    wind_sum.lf=0
    for owpp in buses
        wind_sum.pu=(wind_sum.pu.+deepcopy(owpp.wnd.pu))
        wind_sum.ce=(wind_sum.ce.+deepcopy(owpp.wnd.ce))
        wind_sum.delta=(wind_sum.delta+deepcopy(owpp.wnd.delta))
        wind_sum.lf=(wind_sum.lf+deepcopy(owpp.wnd.lf))
        power_sum=power_sum+deepcopy(owpp.mva)
        push!(paths,deepcopy(as_Astar(domain_oss[owpp.node.num],oss_node,domain_oss)))
    end
    wind_sum.pu=(wind_sum.pu)./length(buses)
    wind_sum.ce=(wind_sum.ce)./length(buses)
    wind_sum.delta=(wind_sum.delta)/length(buses)
    wind_sum.lf=(wind_sum.lf)/length(buses)
    #pick first node in each owpp path
    xys=xy[]
    for path=1:length(paths)
        if (paths[path].parent.num != paths[path].parent.goal)
            push!(xys,deepcopy(paths[path].parent.xy))
        else
            push!(xys,deepcopy(paths[path].xy))
        end
    end

    #re calculate the path
    paths=node[]
    oss_node=opt_OssOptimalLocal(xys,ocn.constrain.ellipses,owp,length(ocn.owpps))
    domain_oss=deepcopy(ocn.discretedom.nodes)
    edges_oss=deepcopy(ocn.discretedom.edges)
    bs,oss_node=opt_ifIn2Out(oss_node,ocn)
    if bs==true
        oss_node.num=deepcopy(length(domain_oss)+1)
        push!(domain_oss,oss_node)
        lof_edgeifyOss(oss_node,ocn,domain_oss,edges_oss)
    end
    return domain_oss,oss_node,power_sum,wind_sum
end

function opt_mvhvOss1stLocal(buses,pcc,ocn)
    #each owpp path
    paths=node[]
    for owp in buses
        push!(paths,deepcopy(as_Astar(owp.node,pcc.node,ocn.discretedom.nodes)))
    end

    #pick first node in each owpp path
    xys=xy[]
    push!(xys,pcc.node.xy)
    for path=1:length(paths)
        pathsParent=paths[path].parent
        while pathsParent.parent.num != buses[path].node.num
            pathsParent=pathsParent.parent
        end
        push!(xys,deepcopy(pathsParent.xy))
    end
    return xys
end

function opt_ttlMvCirc(circ)
    circ.cost=0
    for cb in circ.cbls
        circ.cost=circ.cost+cb.costs.ttl
    end
    for xf in circ.xfmrs
        circ.cost=circ.cost+xf.costs.ttl
    end
end

function opt_ifIn2Out(oss_node,ocn)
    bs,a=lof_test4pnt(oss_node.xy.x,oss_node.xy.y,ocn)
    dummy_dist=Float64
    bsf_dist=Inf
    bsf_indx=1
    if bs==false
        for (indx,pnt) in enumerate(a.nodes)
            dummy_dist=lof_pnt2pnt_dist(oss_node.xy,pnt.xy)
            if dummy_dist<bsf_dist
                bsf_dist=deepcopy(dummy_dist)
                bsf_indx=deepcopy(indx)
            end
        end
        oss_node=a.nodes[bsf_indx]
    else
    end
    return bs,oss_node
end

function opt_buses2nodes4MVoptLocal(mv_square,all,buses)
    target_owpp=1
    owp=bus()
    owp=buses[target_owpp]
    for (indx,bs) in enumerate(buses)
        if buses[target_owpp].node.xy.y<mv_square.ymn
            target_owpp=deepcopy(indx)
            owp=buses[target_owpp]
        end
    end
    return owp
end

function opt_OssOptimalLocal(xys,constrain,owp,num_owpp)
    m = Model(with_optimizer(Ipopt.Optimizer, print_level=0))
    @variable(m, x)
    @variable(m, y)
    @variable(m, lamda)
    eps=1e-9

    #@NLconstraint(m, sum((x-xys[i].x)/(sqrt((xys[i].x-x)^2+(xys[i].y-y)^2)+eps) for i in 1:length(xys))==2*x*lamda)
    #@NLconstraint(m, sum((y-xys[i].y)/(sqrt((xys[i].x-x)^2+(xys[i].y-y)^2)+eps) for i in 1:length(xys))==2*y*lamda)
    @NLobjective(m, Min,sum((sqrt((xys[i].x-x)^2+(xys[i].y-y)^2)+eps) for i in 1:length(xys)))
    #@NLconstraint(m, (x-constrain.x0)^2/(constrain.rx)^2 + (y-constrain.y0)^2/(constrain.ry)^2 == 1)
    #@NLconstraint(m, (x-constrain.x0)^2/(constrain.rx)^2 + (y-constrain.y0)^2/(constrain.ry)^2 >= 1)
    #@NLconstraint(m, (x-constrain.x0)^2/(owpp.zone.pos_width)^2 + (y-constrain.y0)^2/(constrain.ry)^2 <= 1)
    for ellipse in constrain[1:length(constrain)-num_owpp]
        #println(ellipse)
        @NLconstraint(m, (x-ellipse.x0)^2/(ellipse.rx)^2 + (y-ellipse.y0)^2/(ellipse.ry)^2 >= 1)
    end
        @NLconstraint(m, (x-constrain[length(constrain)-num_owpp+owp.num].x0)^2/(constrain[length(constrain)-num_owpp+owp.num].rx)^2 + (y-constrain[length(constrain)-num_owpp+owp.num].y0)^2/(constrain[length(constrain)-num_owpp+owp.num].ry)^2 >= 1)
    #println(mv_square.ymx)
    #=println(mv_square.ymn)
    println(mv_square.xmx)
    println(mv_square.xmn)

    println(xys)=#
    #@constraint(m, x >= mv_square.xmn)
    @constraint(m, y == owp.node.xy.y)
    #@constraint(m, x <= mv_square.xmx)
    #@constraint(m, y <= mv_square.ymx)
    optimize!(m)
    temp_xy=xy()
    temp_xy.x=JuMP.value.((x))
    temp_xy.y=JuMP.value.((y))
    node_oss=node()
    node_oss.xy=temp_xy

    #node_oss=opt_findClosestNode(temp_xy,owpp.zone.nodes)
    return node_oss
end

function opt_findClosestNode(x_y,ndes)
    node_oss=node()
    distance=Inf
    oss_indx=Inf

    for (indx,nd) in enumerate(ndes)
        temp_dist=lof_pnt2pnt_dist(x_y,nd.xy)
        if(temp_dist<distance)
            distance=deepcopy(temp_dist)
            oss_indx=deepcopy(indx)
        end
    end
    node_oss=ndes[oss_indx]
    return node_oss
end

#MV constraints
function opt_mvConstraints(ocn, owpps)
    sqr_bnd=square()
    sqr_bnd.ymn=0
    sqr_bnd.ymx=Inf
    sqr_bnd.xmn=0
    sqr_bnd.xmx=Inf
    for (indx,owpp) in enumerate(ocn.owpps)

        #Check MV ranges
        if (owpps[indx]==1)
            ymn=owpp.node.xy.y-owpp.mv_zone.neg_height
            ymx=owpp.node.xy.y+owpp.mv_zone.pos_height
            xmn=owpp.node.xy.x-owpp.mv_zone.neg_width
            xmx=owpp.node.xy.x+owpp.mv_zone.pos_width

            #south boundary
            if (ymn>sqr_bnd.ymn)
                sqr_bnd.ymn=deepcopy(ymn)
            end

            #north boundary
            if (ymx<sqr_bnd.ymx)
                sqr_bnd.ymx=deepcopy(ymx)
            end

            #west boundary
            if (xmn>sqr_bnd.xmn)
                sqr_bnd.xmn=deepcopy(xmn)
            end

            #east boundary
            if (xmx<sqr_bnd.xmx)
                sqr_bnd.xmx=deepcopy(xmx)
            end
        else
        end
    end
    return sqr_bnd
end

#OWWP constraints
function opt_owppConstraints(ocn)
    #@NLconstraint(m, ((x-x1)/a)^2+((y-y1)/b)^2 > 1)
    for owpp in ocn.owpps
        ellipse_con=ellipse()
        ellipse_con.y0=owpp.node.xy.y-(owpp.zone.neg_height-owpp.zone.pos_height)/2
        ellipse_con.x0=owpp.node.xy.x-(owpp.zone.neg_width-owpp.zone.pos_width)/2
        ellipse_con.ry=owpp.node.xy.y+owpp.zone.pos_height-ellipse_con.y0
        ellipse_con.rx=owpp.node.xy.x+owpp.zone.pos_width-ellipse_con.x0
        push!(ocn.constrain.ellipses,deepcopy(ellipse_con))
    end
end


#OWWP constraints
function opt_nogoConstraints(ocn)
    #@NLconstraint(m, ((x-x1)/a)^2+((y-y1)/b)^2 > 1)
    for nogo in ocn.nogos
        ellipse_con=ellipse()
        cntr=opt_findShapeCentre(nogo)
        ellipse_con.y0=cntr.y
        ellipse_con.x0=cntr.x
        ellipse_con.ry=opt_findShapeRadius(nogo,cntr)
        ellipse_con.rx=ellipse_con.ry
        push!(ocn.constrain.ellipses,deepcopy(ellipse_con))
    end
end
#=
ng=ocean.nogos[1]
=#
#Finds closest edge
function opt_findShapeRadius(ng,cntr)
    rad=Inf
    for bnds in [ng.nbnd,ng.sbnd,ng.ebnd,ng.wbnd]
        for bnd in bnds
            x=(cntr.x+(cntr.y*bnd.m_findx+bnd.b_findx))/2
            y=(cntr.y+(cntr.x*bnd.m_findy+bnd.b_findy))/2
            rad_dummy=sqrt((cntr.x-x)^2+(cntr.y-y)^2)
            if (rad_dummy<rad)
                rad=deepcopy(rad_dummy)
            end
        end
    end
    return rad
end

#=
ocean.owpps[2].node
mean=opt_findShapeCentre(ocean.owpps[2].zone)
=#
function opt_findShapeCentre(shape_bnds)
    mn0=opt_shapeMean(shape_bnds)
    triangles=Array{xy,1}()
    for bnd in shape_bnds.nbnd
        xy0=xy()
        xy1=xy()
        xy0.x=bnd.xmn
        xy0.y=bnd.xmn*bnd.m_findy+bnd.b_findy
        xy1.x=bnd.xmx
        xy1.y=bnd.xmx*bnd.m_findy+bnd.b_findy
        tri_mean=opt_findTriAngCentre(mn0, xy0, xy1)
        push!(triangles,deepcopy(tri_mean))
    end

    for bnd in shape_bnds.sbnd
        xy0=xy()
        xy1=xy()
        xy0.x=bnd.xmn
        xy0.y=bnd.xmn*bnd.m_findy+bnd.b_findy
        xy1.x=bnd.xmx
        xy1.y=bnd.xmx*bnd.m_findy+bnd.b_findy
        tri_mean=opt_findTriAngCentre(mn0, xy0, xy1)
        push!(triangles,deepcopy(tri_mean))
    end

    for bnd in shape_bnds.ebnd
        xy0=xy()
        xy1=xy()
        xy0.x=bnd.ymn*bnd.m_findx+bnd.b_findx
        xy0.y=bnd.ymn
        xy1.x=bnd.ymx*bnd.m_findx+bnd.b_findx
        xy1.y=bnd.ymx
        tri_mean=opt_findTriAngCentre(mn0, xy0, xy1)
        push!(triangles,deepcopy(tri_mean))
    end

    for bnd in shape_bnds.wbnd
        xy0=xy()
        xy1=xy()
        xy0.x=bnd.ymn*bnd.m_findx+bnd.b_findx
        xy0.y=bnd.ymn
        xy1.x=bnd.ymx*bnd.m_findx+bnd.b_findx
        xy1.y=bnd.ymx
        tri_mean=opt_findTriAngCentre(mn0, xy0, xy1)
        push!(triangles,deepcopy(tri_mean))
    end

    cntr=xy()
    cntr.x=0
    cntr.y=0
    for x_y in triangles
        cntr.x=cntr.x+deepcopy(x_y.x)
        cntr.y=cntr.y+deepcopy(x_y.y)
    end
    cntr.x=cntr.x/length(triangles)
    cntr.y=cntr.y/length(triangles)
    return cntr
end

function opt_shapeMean(shape_bnds)
    mean=xy()
    mean.y=0
    mean.x=0

    #north bnd
    for bnd in shape_bnds.nbnd
        mean.y=mean.y+deepcopy(bnd.ymn)+deepcopy(bnd.ymx)
        mean.x=mean.x+deepcopy(bnd.xmn)+deepcopy(bnd.xmx)
    end

    #south bnd
    for bnd in shape_bnds.sbnd
        mean.y=mean.y+deepcopy(bnd.ymn)+deepcopy(bnd.ymx)
        mean.x=mean.x+deepcopy(bnd.xmn)+deepcopy(bnd.xmx)
    end
    #west bnd
    for bnd in shape_bnds.wbnd
        mean.y=mean.y+deepcopy(bnd.ymn)+deepcopy(bnd.ymx)
        mean.x=mean.x+deepcopy(bnd.xmn)+deepcopy(bnd.xmx)
    end
    #east bnd
    for bnd in shape_bnds.ebnd
        mean.y=mean.y+deepcopy(bnd.ymn)+deepcopy(bnd.ymx)
        mean.x=mean.x+deepcopy(bnd.xmn)+deepcopy(bnd.xmx)
    end
    l=2*(length(shape_bnds.ebnd)+length(shape_bnds.wbnd)+length(shape_bnds.sbnd)+length(shape_bnds.nbnd))
    mean.y=mean.y/(l)
    mean.x=mean.x/(l)
    return mean
end
#=
y0=ocean.nogos[1].wbnd[1].ymn
y1=ocean.nogos[1].wbnd[1].ymx
y2=ocean.nogos[1].ebnd[1].ymn
x0=ocean.nogos[1].wbnd[1].m_findx*y0+ocean.nogos[1].wbnd[1].b_findx
x1=ocean.nogos[1].wbnd[1].m_findx*y1+ocean.nogos[1].wbnd[1].b_findx
x2=ocean.nogos[1].ebnd[1].m_findx*y2+ocean.nogos[1].ebnd[1].b_findx
pnt0=xy()
pnt1=xy()
pnt2=xy()
pnt0.x=x0
pnt0.y=y0
pnt1.x=x1
pnt1.y=y1
pnt2.x=x2
pnt2.y=y2
cntr=opt_findTriAngCentre(pnt0, pnt1, pnt2)
=#
function opt_findTriAngCentre(p0, p1, p2)
    cntr=xy()
    cntr.x=(p0.x+p1.x+p2.x)/3
    cntr.y=(p0.y+p1.y+p2.y)/3
    return cntr
end
