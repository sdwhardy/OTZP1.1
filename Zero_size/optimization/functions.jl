include("compound/functions.jl")#

################################### MV Topologies ##############################
################################################################################

#MV topology main function
function opt_mvOSSplacement(ocn,pcc)
    owpps_tbl=top_hvTopos(length(ocn.owpps))
    owpps_tbl=owpps_tbl[end:-1:1,end:-1:1]
    bus_dummies=Array{bus,1}()
    oss_systems=Array{Array{circuit, 1}, 1}()
    dummy_sys=circuit()
    dummy_sys.cost=0

    for indx0=1:length(owpps_tbl[:,1])
        bus_dummies=bus[]
        validMV,mv_regionOWPP=opt_mvConstraints(ocn,owpps_tbl[indx0,:])

        if (validMV==true)
            for (indx1,owpp) in enumerate(owpps_tbl[indx0,:])
                if (owpp==1)
                    push!(bus_dummies,ocn.owpps[indx1])
                end
            end
            oss_system=opt_mvOssSystem(mv_regionOWPP,bus_dummies,pcc,ocn,owpps_tbl[indx0,:])
            oss_system.id=string(oss_system.decimal)*"#m"
            push!(oss_systems,[oss_system])
        else
            push!(oss_systems,[dummy_sys])
        end
    end
    return oss_systems
end

function opt_mvConstraints(ocn, owpps)
    mv_regionOWPP=bus[]
    valid=false
    if (sum(owpps)>1)
        first_owp=findfirst(x->x==1,owpps)
        push!(mv_regionOWPP,ocn.owpps[first_owp])
        for (i,owp) in enumerate(owpps[first_owp+1:length(owpps)])
            if (owp==1)
                mv_rng=true
                for mv_owp in mv_regionOWPP
                    o2oL=lof_pnt2pnt_dist(mv_owp.node.xy,ocn.owpps[i+first_owp].node.xy)
                    if (o2oL>=mv_owp.mv_zone+ocn.owpps[i+first_owp].mv_zone)
                        mv_rng=false
                    end
                end
                if (mv_rng==true)
                    push!(mv_regionOWPP,ocn.owpps[i+first_owp])
                end
            end
        end
    end
    if (length(mv_regionOWPP)>1)
        valid=true
    else
        valid=false
    end
    return valid,mv_regionOWPP
end
#buses=bus_dummies
#bn=owpps_tbl[indx0,:]
function opt_mvOssSystem(mv_regionOWPP,buses,pcc,ocn,bn)
    circ=circuit()
    opNode = opt_findOptPointMV(buses,ocn,mv_regionOWPP)
    circ=opt_mvhvEquipment(circ,buses[1],buses,pcc,ocn,bn,opNode)
    return circ
end

function opt_findOptPointMV(buses,ocn,mv_regionOWPP)
    rmnXys=Array{xy,1}()
    for bs in buses
        push!(rmnXys,deepcopy(bs.node.xy))
    end
    opNode=opt_OssOptimalMV(rmnXys,buses,mv_regionOWPP)
    ocn.buses=ocn.buses+1
    opNode.num=deepcopy(ocn.buses)
    return opNode
end

function opt_OssOptimalMV(xys,buses,mv_regionOWPP)
    m = Model(with_optimizer(Ipopt.Optimizer, print_level=0))
    @variable(m, x)
    @variable(m, y)
    eps=1e-8

    @NLobjective(m, Min,sum(sqrt((xys[i].x-x)^2+(xys[i].y-y)^2+eps) for i in 1:length(xys)))
    for (ind,op) in enumerate(mv_regionOWPP)
        @NLconstraint(m, (x-op.node.xy.x)^2+(y-op.node.xy.y)^2 <= ((op.mv_zone-(eps))^2))
    end

    optimize!(m)
    temp_xy=xy()
    temp_xy.x=JuMP.value.((x))
    temp_xy.y=JuMP.value.((y))
    node_oss=node()
    node_oss.xy=temp_xy
    return node_oss
end

function opt_updateMVocean(ocn)
    for (i,circ) in enumerate(ocn.mv_circuits)
        if (circ[length(circ)].cost<1)
            ocn.mv_circuits[i]=ocn.hv_circuits[i]
        end
    end
    return ocn.mv_circuits
end
################################### HV Topologies ##############################
################################################################################
#**
function opt_hvOSSplacement(ocn,pcc)
    hv_circs=Array{Array{circuit, 1}, 1}()
    owpps_tbl_hv=top_hvTopos(length(ocn.owpps))
    owpps_tbl_hv=owpps_tbl_hv[end:-1:1,end:-1:1]
    oss_system=circuit()
    oss_systems=Array{circuit,1}()
    for indx0=1:length(owpps_tbl_hv[:,1])
        bus_dummies=bus[]
        for (indx1,owp) in enumerate(owpps_tbl_hv[indx0,:])
            if (owp==1)
                push!(bus_dummies,ocn.owpps[indx1])
            end
        end

        if (length(bus_dummies)>1)
            oss_system=opt_hvOssSystem(bus_dummies,pcc,ocn,owpps_tbl_hv[indx0,:])
        else
            oss_system=opt_str8Connect(bus_dummies[1],pcc,ocn,owpps_tbl_hv[indx0,:])
        end
        oss_system.id=string(oss_system.decimal)*"#h"
        push!(hv_circs,[oss_system])
    end
    return hv_circs
end

#single owpp to pcc connection, used in HV connection
function opt_str8Connect(owp,pcc,ocn,bn)
    circ=circuit()
    circ.binary=bn
    circ.decimal=top_bin2dec(bn)
    circ.pcc=deepcopy(pcc)
    circ.base_owp=owp
    push!(circ.owpps,owp)

    #Calculate length and type of cables
    cs=cstF_MvHvCblpcc(lof_pnt2pnt_dist(owp.node.xy,pcc.node.xy),owp.mva,owp.wnd,ocn.finance,pcc,ocn.sys,owp.mv_zone,ocn.eqp_data)
    #store cable data
    if length(cs) == 2
        push!(cs[1].pth,deepcopy(owp.node))
        push!(cs[1].pth,deepcopy(pcc.node))
        push!(circ.pcc_cbls,deepcopy(cs[1]))
        #PCC
        push!(circ.pcc.xfmrs,deepcopy(cs[2]))
    elseif length(cs) == 4
        mvhvOSSnd=opt_oss4mv2hv(pcc.node,owp)
        ocn.buses=ocn.buses+1
        mvhvOSSnd.num=ocn.buses
        #Cable
        push!(cs[1].pth,deepcopy(owp.node))
        push!(cs[1].pth,deepcopy(mvhvOSSnd))
        push!(cs[2].pth,deepcopy(mvhvOSSnd))
        push!(cs[2].pth,deepcopy(pcc.node))

        push!(circ.owp_MVcbls,deepcopy(cs[1]))
        push!(circ.pcc_cbls,deepcopy(cs[2]))

        #OSS
        ossmv=bus()
        ossmv.node=mvhvOSSnd
        ossmv.base_cost=ocn.finance.FC_bld
        ossmv.mva=deepcopy(owp.mva)
        ossmv.wnd=deepcopy(owp.wnd)
        push!(circ.osss_mog,ossmv)
        push!(circ.osss_mog[1].xfmrs,deepcopy(cs[3]))
        #PCC
        push!(circ.pcc.xfmrs,deepcopy(cs[4]))
    end
    circ.oss_mva=deepcopy(owp.mva)
    circ.oss_wind=deepcopy(owp.wnd)
    opt_ttlMvCirc(circ)
    return circ
end


function opt_hvOssSystem(buses,pcc,ocn,bn)
    circ=circuit()
    owp=buses[1]
    opNode = opt_findOptPointHV(buses,pcc,ocn,owp)
    circ=opt_mvhvEquipment(circ,owp,buses,pcc,ocn,bn,opNode)
    return circ
end


#find the optimal position of an HV connected topology OSS
function opt_findOptPointHV(buses,pcc,ocn,owp)
    rmnxys=Array{xy,1}()
    for bs in buses
        push!(rmnxys,deepcopy(bs.node.xy))
    end
    opNode=opt_HVOssOptimal(rmnxys,owp,pcc)
    ocn.buses=ocn.buses+1
    opNode.num=deepcopy(ocn.buses)
    return opNode
end

#xys=rmnxys
function opt_HVOssOptimal(xys,owp,pcc)
    m = Model(with_optimizer(Ipopt.Optimizer, print_level=1))
    @variable(m, x)
    @variable(m, y)
    x_mn=pcc.node.xy.x    #x_mn=18.215405
    x_mx=pcc.node.xy.x#x_mx=22.73574
    y_mn=pcc.node.xy.y
    y_mx=owp.node.xy.y
    for co_ord in xys
        if (co_ord.x<x_mn)
            x_mn=deepcopy(co_ord.x)
        elseif (co_ord.x>x_mx)
            x_mx=deepcopy(co_ord.x)
        end
        if (co_ord.y<y_mx)
            y_mx=deepcopy(co_ord.y)
        end
    end
    #@variable(m, lamda)
    epsil=1e-9
    push!(xys,pcc.node.xy)
    @NLobjective(m, Min,sum(sqrt((xys[i].x-x)^2+(xys[i].y-y)^2+epsil) for i in 1:length(xys)))

    @constraint(m, y <= y_mx)
    @constraint(m, y >= y_mn)
    @constraint(m, x <= x_mx)
    @constraint(m, x >= x_mn)

    optimize!(m)
    temp_xy=xy()
    temp_xy.x=JuMP.value.((x))
    temp_xy.y=JuMP.value.((y))
    node_oss=node()
    node_oss.xy=temp_xy
    return node_oss
end


#**
#List all combinatorial connections of HV topologies (returns a binary table)
function top_hvTopos(Numba)
    #find all combinations
    clms=trunc(Int,Numba)
    rows=trunc(Int, 2.0^clms)
    empty_tbl=eensF_blankTbl(rows,clms)

    #remove anything less than a pair
    lngth=length(empty_tbl[:,1])
    indx=1
    while indx <= lngth
        if (sum(empty_tbl[indx,:]) < 1)
            empty_tbl=empty_tbl[1:size(empty_tbl,1) .!= indx,: ]
            indx=indx-1
        end
        indx=indx+1
        lngth=length(empty_tbl[:,1])
    end
    return empty_tbl
end
######################### Optimize topology ##########################
######################################################################

function opt_readjust_circuits(ocn,circs)
    longest_cable=lof_pnt2pnt_dist(ocn.owpps[length(ocn.owpps)].node.xy,ocn.pccs[length(ocn.pccs)].node.xy)
    for (i0,cs0) in enumerate(circs)
        for (i1,cs1) in enumerate(cs0)
            c_o=deepcopy(circs[i0][i1])
            new_coords=opt_reAdjust_oss(circs[i0][i1],ocn.owpps[1].mv_zone,ocn.sys.mvCl,10e-6)
            circs[i0][i1]=opt_reAdjust_cbls(circs[i0][i1],new_coords,ocn,longest_cable)
            if ((c_o.cost)<circs[i0][i1].cost)
                #println("Error: re-adjustment failed, circuit: ")
                #print(string(i)*" Initial: "*string(c_o))
                #print(" - Adjusted: "*string(circs[i].cost))
                new_coords=opt_reAdjust_oss(circs[i0][i1],ocn.owpps[1].mv_zone,ocn.sys.mvCl,10e-8)
                circs[i0][i1]=opt_reAdjust_cbls(circs[i0][i1],new_coords,ocn,longest_cable)
                #print(" - re-adjusted: "*string(circs[i].cost))
                if ((c_o.cost)<circs[i0][i1].cost)
                    circs[i0][i1]=deepcopy(c_o)
                    #println("kept original layout!")
                end
            end
        end
    end
    return circs
end


function opt_reAdjust_oss(system,mv_rng,mv_cl,eps)
    mog_xys=Array{Tuple{xy,Int64},1}()
    oss_xys=Array{Tuple{xy,Int64},1}()
    for mog in system.osss_mog
        push!(mog_xys,(mog.node.xy,mog.node.num))
    end
    for os in system.osss_owp
        push!(oss_xys,(os.node.xy,os.node.num))
    end
    m = Model(with_optimizer(Ipopt.Optimizer, print_level=1))

    @variables(m, begin
               MOG_y[i=1:length(mog_xys)]>=0, (start = mog_xys[i][1].y, base_name = "MOGY_"*string(mog_xys[i][2]))
               MOG_x[i=1:length(mog_xys)]>=0, (start = mog_xys[i][1].x, base_name = "MOGX_"*string(mog_xys[i][2]))
           end)

       @variables(m, begin
                  OSS_y[i=1:length(oss_xys)]>=0, (start = oss_xys[i][1].y, base_name = "OSSY_"*string(oss_xys[i][2]))
                  OSS_x[i=1:length(oss_xys)]>=0, (start = oss_xys[i][1].x, base_name = "OSSX_"*string(oss_xys[i][2]))
              end)

    connections=Array{Tuple{Any,Any,Any,Any,Float64},1}()
    for mv_c in system.owp_MVcbls
        for (i,mg) in enumerate(mog_xys)
            if (mg[2]==mv_c.pth[length(mv_c.pth)].num)
                push!(connections,(mv_c.pth[1].xy.x, mv_c.pth[1].xy.y, MOG_x[i], MOG_y[i], mv_c.costs.ttl/mv_c.length))
                @constraint(m,((connections[length(connections)][1]-connections[length(connections)][3])^2+(connections[length(connections)][2]-connections[length(connections)][4])^2) <= mv_rng^2)
            end
        end
        for (i,os) in enumerate(oss_xys)
            if (os[2]==mv_c.pth[length(mv_c.pth)].num)
                push!(connections,(mv_c.pth[1].xy.x, mv_c.pth[1].xy.y, OSS_x[i], OSS_y[i], mv_c.costs.ttl/mv_c.length))
                @constraint(m,((connections[length(connections)][1]-connections[length(connections)][3])^2+(connections[length(connections)][2]-connections[length(connections)][4])^2) == mv_cl^2)
            end
        end
    end

    #PCC_connections=Array{Tuple{VariableRef,VariableRef,Float64,Float64,Float64},1}()
    for pcc_c in system.pcc_cbls
        for (i,mg) in enumerate(mog_xys)
            if (mg[2]==pcc_c.pth[1].num)
                push!(connections,(MOG_x[i], MOG_y[i], pcc_c.pth[length(pcc_c.pth)].xy.x, pcc_c.pth[length(pcc_c.pth)].xy.y, pcc_c.costs.ttl/pcc_c.length))
            end
        end
        for (i,os) in enumerate(oss_xys)
            if (os[2]==pcc_c.pth[length(pcc_c.pth)].num)
                push!(connections,(OSS_x[i], OSS_y[i], pcc_c.pth[length(pcc_c.pth)].xy.x, pcc_c.pth[length(pcc_c.pth)].xy.y, pcc_c.costs.ttl/pcc_c.length))
            end
        end
    end

    #HV_connections=Array{Tuple{VariableRef,VariableRef,VariableRef,VariableRef,Float64},1}()
    for hv_c in system.owp_HVcbls
        tail_x=VariableRef
        tail_y=VariableRef
        head_x=VariableRef
        head_y=VariableRef
        for (i,mg) in enumerate(mog_xys)

            if (mg[2]==hv_c.pth[length(hv_c.pth)].num)
                head_x=MOG_x[i]
                head_y=MOG_y[i]
            end
        end
        for (i,os) in enumerate(oss_xys)
            if (os[2]==hv_c.pth[1].num)
                tail_x=OSS_x[i]
                tail_y=OSS_y[i]
            end
        end
        push!(connections,(tail_x, tail_y, head_x, head_y, hv_c.costs.ttl/hv_c.length))
    end
    #o2o_connections=Array{Tuple{VariableRef,VariableRef,VariableRef,VariableRef,Float64},1}()
    for o2o_c in system.oss2oss_cbls
        tail_x=VariableRef
        tail_y=VariableRef
        head_x=VariableRef
        head_y=VariableRef
        for (i,mg) in enumerate(mog_xys)
            if (mg[2]==o2o_c.pth[1].num)
                tail_x=MOG_x[i]
                tail_y=MOG_y[i]
            end
            if (mg[2]==o2o_c.pth[length(o2o_c.pth)].num)
                head_x=MOG_x[i]
                head_y=MOG_y[i]
            end
        end

        push!(connections,(tail_x, tail_y, head_x, head_y, o2o_c.costs.ttl/o2o_c.length))
    end

    @NLobjective(m, Min,sum(sqrt((((connections[i][1]-connections[i][3])^2+(connections[i][2]-connections[i][4])^2)+eps))*connections[i][5] for i in 1:length(connections)))

    optimize!(m)
    mxy=(JuMP.value.(MOG_x),JuMP.value.(MOG_y),JuMP.value.(OSS_x),JuMP.value.(OSS_y))
    return mxy
end
#circ=circs[i]
#co_ords=new_coords
function opt_reAdjust_cbls(circ,co_ords,ocn,longest_cable)
    for mg_i=1:1:length(circ.osss_mog)
        circ.osss_mog[mg_i].node.xy.x=deepcopy(co_ords[1][mg_i])
        circ.osss_mog[mg_i].node.xy.y=deepcopy(co_ords[2][mg_i])
    end
    for oss_i=1:1:length(circ.osss_owp)
        circ.osss_owp[oss_i].node.xy.x=deepcopy(co_ords[3][oss_i])
        circ.osss_owp[oss_i].node.xy.y=deepcopy(co_ords[4][oss_i])
    end
    #set
    circ.owp_MVcbls=opt_updateMVC(circ,ocn)
    circ.owp_HVcbls=opt_updateHVC(circ,ocn,longest_cable)
    circ.oss2oss_cbls=opt_updateo2o(circ,ocn,longest_cable)
    circ.pcc_cbls=opt_updatePcc(circ,ocn,2*longest_cable)
    opt_ttlMvCirc(circ)
    return circ
end

function opt_updateMVC(circ,ocn)
    for mvc_i=1:1:length(circ.owp_MVcbls)
        mvc_owpNodeNum=deepcopy(circ.owp_MVcbls[mvc_i].pth[1].num)
        mvc_ossNodeNum=deepcopy(circ.owp_MVcbls[mvc_i].pth[length(circ.owp_MVcbls[mvc_i].pth)].num)
        circ.owp_MVcbls[mvc_i].pth=node[]
        tale=bus()
        hed=bus()
        kble=cbl()
        for opp in circ.owpps
            if (opp.node.num==mvc_owpNodeNum)
                tale=deepcopy(opp)
                break
            end
        end
        for oss in circ.osss_owp
            if (oss.node.num==mvc_ossNodeNum)
                hed=deepcopy(oss)
                break
            end
        end
        for mog in circ.osss_mog
            if (mog.node.num==mvc_ossNodeNum)
                hed=deepcopy(mog)
                break
            end
        end
        mv_l=lof_pnt2pnt_dist(tale.node.xy,hed.node.xy)
        if (mv_l<circ.owp_MVcbls[mvc_i].length)
            if (mv_l<ocn.sys.mvCl)
                mv_l=ocn.sys.mvCl
            end
            kble=cstF_MvCbl_nextSizeDown(mv_l,tale.mva,circ.owp_MVcbls[mvc_i].elec.volt,tale.wnd,ocn.finance,circ.owp_MVcbls[mvc_i].size,circ.owp_MVcbls[mvc_i].num,ocn.eqp_data)
            circ.owp_MVcbls[mvc_i]=kble
        elseif (mv_l>circ.owp_MVcbls[mvc_i].length)
            if (mv_l<ocn.owpps[length(ocn.owpps)].mv_zone+10e-3)
                kble=cstF_MvCbl_nextSizeUp(mv_l,tale.mva,circ.owp_MVcbls[mvc_i].elec.volt,tale.wnd,ocn.finance,circ.owp_MVcbls[mvc_i].size,circ.owp_MVcbls[mvc_i].num,ocn.eqp_data)
                circ.owp_MVcbls[mvc_i]=kble
            else
                circ.owp_MVcbls[mvc_i].costs.ttl=Inf
            end
        else
        end
        push!(circ.owp_MVcbls[mvc_i].pth,tale.node)
        push!(circ.owp_MVcbls[mvc_i].pth,hed.node)
    end
    return circ.owp_MVcbls
end

function opt_updateHVC(circ,ocn,longest_cable)
    S=0
    wp=wind()
    for hvc_i=1:1:length(circ.owp_HVcbls)
        hvc_ossNodeNum=deepcopy(circ.owp_HVcbls[hvc_i].pth[1].num)
        hvc_mogNodeNum=deepcopy(circ.owp_HVcbls[hvc_i].pth[length(circ.owp_HVcbls[hvc_i].pth)].num)
        circ.owp_HVcbls[hvc_i].pth=node[]
        tale=node()
        hed=node()
        kble=cbl()
        #println(hvc_ossNodeNum)
        for oss in circ.osss_owp
            #println(oss.node.goal)
            if (oss.node.num==hvc_ossNodeNum)
                tale=oss.node
                S=oss.mva
                wp=oss.wnd
                break
            end
        end
        for mog in circ.osss_mog
            if (mog.node.num==hvc_mogNodeNum)
                hed=mog.node
                break
            end
        end
        hv_l=lof_pnt2pnt_dist(tale.xy,hed.xy)
        if (hv_l<circ.owp_HVcbls[hvc_i].length)
            kble=cstF_nextSizeDown(hv_l,S,circ.owp_HVcbls[hvc_i].elec.volt,wp,ocn.finance,circ.owp_HVcbls[hvc_i].size,circ.owp_HVcbls[hvc_i].num,ocn.eqp_data)
            circ.owp_HVcbls[hvc_i]=kble
        elseif (hv_l>circ.owp_HVcbls[hvc_i].length)
            if (hv_l<longest_cable)
                kble=cstF_nextSizeUp(hv_l,S,circ.owp_HVcbls[hvc_i].elec.volt,wp,ocn.finance,circ.owp_HVcbls[hvc_i].size,circ.owp_HVcbls[hvc_i].num,ocn.eqp_data)
                circ.owp_HVcbls[hvc_i]=kble
            else
                circ.owp_HVcbls[hvc_i].costs.ttl=Inf
            end
        else
        end
        push!(circ.owp_HVcbls[hvc_i].pth,tale)
        push!(circ.owp_HVcbls[hvc_i].pth,hed)
    end
    return circ.owp_HVcbls
end

function opt_updateo2o(circ,ocn,longest_cable)
    S=0
    wp=wind()
    for o2o_i=1:1:length(circ.oss2oss_cbls)
        o2o_ossNodeNum=deepcopy(circ.oss2oss_cbls[o2o_i].pth[1].num)
        o2o_mogNodeNum=deepcopy(circ.oss2oss_cbls[o2o_i].pth[length(circ.oss2oss_cbls[o2o_i].pth)].num)
        circ.oss2oss_cbls[o2o_i].pth=node[]
        tale=node()
        hed=node()
        kble=cbl()
        for mog in circ.osss_mog
            if (mog.node.num==o2o_ossNodeNum)
                tale=deepcopy(mog.node)
                S=deepcopy(mog.mva)
                wp=deepcopy(mog.wnd)
            elseif (mog.node.num==o2o_mogNodeNum)
                hed=deepcopy(mog.node)
            else
            end
        end
        o2o_l=lof_pnt2pnt_dist(tale.xy,hed.xy)
        if (o2o_l<circ.oss2oss_cbls[o2o_i].length)
            kble=cstF_nextSizeDown(o2o_l,S,circ.oss2oss_cbls[o2o_i].elec.volt,wp,ocn.finance,circ.oss2oss_cbls[o2o_i].size,circ.oss2oss_cbls[o2o_i].num,ocn.eqp_data)
            circ.oss2oss_cbls[o2o_i]=kble
        elseif (o2o_l>circ.oss2oss_cbls[o2o_i].length)
            if (o2o_l<longest_cable)
                kble=cstF_nextSizeUp(o2o_l,S,circ.oss2oss_cbls[o2o_i].elec.volt,wp,ocn.finance,circ.oss2oss_cbls[o2o_i].size,circ.oss2oss_cbls[o2o_i].num,ocn.eqp_data)
                circ.oss2oss_cbls[o2o_i]=kble
            else
                circ.oss2oss_cbls[o2o_i].costs.ttl=Inf
            end
        else
        end
        push!(circ.oss2oss_cbls[o2o_i].pth,tale)
        push!(circ.oss2oss_cbls[o2o_i].pth,hed)
    end
    return circ.oss2oss_cbls
end

function opt_updatePcc(circ,ocn,longest_cable)
    S=0
    wp=wind()
    for pcc_i=1:1:length(circ.pcc_cbls)
        pcc_mogNodeNum=deepcopy(circ.pcc_cbls[pcc_i].pth[1].num)
        circ.pcc_cbls[pcc_i].pth=node[]
        tale=node()
        hed=ocn.pccs[length(ocn.pccs)].node
        kble=cbl()
        for mog in circ.osss_mog
            if (mog.node.num==pcc_mogNodeNum)
                tale=deepcopy(mog.node)
                S=deepcopy(mog.mva)
                wp=deepcopy(mog.wnd)
            else
            end
        end
        for oss in circ.osss_owp
            if (oss.node.num==pcc_mogNodeNum)
                tale=deepcopy(oss.node)
                S=deepcopy(oss.mva)
                wp=deepcopy(oss.wnd)
            else
            end
        end
        pcc_l=lof_pnt2pnt_dist(tale.xy,hed.xy)
        if (pcc_l<circ.pcc_cbls[pcc_i].length)
            kble=cstF_nextSizeDownPcc(pcc_l,S,circ.pcc_cbls[pcc_i].elec.volt,wp,ocn.finance,circ.pcc_cbls[pcc_i].size,circ.pcc_cbls[pcc_i].num,ocn.eqp_data)
            circ.pcc_cbls[pcc_i]=kble
        elseif (pcc_l>circ.pcc_cbls[pcc_i].length)
            if (pcc_l<longest_cable)
                kble=cstF_nextSizeUpPcc(pcc_l,S,circ.pcc_cbls[pcc_i].elec.volt,wp,ocn.finance,circ.pcc_cbls[pcc_i].size,circ.pcc_cbls[pcc_i].num,ocn.eqp_data)
                circ.pcc_cbls[pcc_i]=kble
            else
                circ.pcc_cbls[pcc_i].costs.ttl=Inf
            end
        else
        end
        push!(circ.pcc_cbls[pcc_i].pth,tale)
        push!(circ.pcc_cbls[pcc_i].pth,hed)
    end
    return circ.pcc_cbls
end

########################### General MV/HV ############################
######################################################################
#Main cost fuction for both MV and HV connections

function opt_mvhvEquipment(circ,owp,buses,pcc,ocn,bn,opNode)
    circ.base_owp=owp
    circ.oss_mva,circ.oss_wind=opt_findWindPwr(buses)

    #Find PCC connection cable and PCC transformer
    cbl_xfo=cstF_HvCblallKvo2p(lof_pnt2pnt_dist(pcc.node.xy,opNode.xy),circ.oss_mva,circ.oss_wind,ocn.finance,pcc,ocn.eqp_data)
    #Set Pcc cable path
    push!(cbl_xfo[1].pth,deepcopy(opNode))
    push!(cbl_xfo[1].pth,deepcopy(pcc.node))
    push!(circ.pcc_cbls,cbl_xfo[1])

    #pcc transformer
    circ.pcc=deepcopy(pcc)
    push!(circ.pcc.xfmrs,cbl_xfo[2])

    #MV cable sizes and paths
    pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400=opt_InitPW()
    for (i,op) in enumerate(buses)
        cb,xs=cstF_MvHvCbloss(lof_pnt2pnt_dist(opNode.xy,op.node.xy),op.mva,op.wnd,ocn.finance,circ.pcc_cbls[1].elec.volt,ocn.sys,op.mv_zone,ocn.eqp_data)
        if (length(cb)==1)
            push!(cb[1].pth,deepcopy(op.node))
            push!(cb[1].pth,deepcopy(opNode))
            push!(circ.owp_MVcbls,cb[1])
            #store winds and powers
            push!(pMv,op.mva)
            push!(wMv,op.wnd)
        elseif (length(cb)==2)
            #Cable
            mvhvOSSnd=opt_oss4mv2hv(opNode,op)
            ocn.buses=ocn.buses+1
            mvhvOSSnd.num=ocn.buses

            push!(cb[1].pth,deepcopy(op.node))
            push!(cb[1].pth,deepcopy(mvhvOSSnd))
            push!(cb[2].pth,deepcopy(mvhvOSSnd))
            push!(cb[2].pth,deepcopy(opNode))
            push!(circ.owp_MVcbls,deepcopy(cb[1]))
            push!(circ.owp_HVcbls,deepcopy(cb[2]))

            #OSS
            ossmv=bus()
            ossmv.node=mvhvOSSnd
            ossmv.base_cost=ocn.finance.FC_bld
            ossmv.mva=op.mva
            ossmv.wnd=op.wnd
            push!(circ.osss_owp,ossmv)
            push!(circ.osss_owp[length(circ.osss_owp)].xfmrs,deepcopy(xs[1]))

            #store winds and powers
            if (circ.owp_HVcbls[length(circ.owp_HVcbls)].elec.volt == circ.pcc_cbls[length(circ.pcc_cbls)].elec.volt)
                push!(pHv,op.mva)
                push!(wHv,op.wnd)
            elseif (circ.owp_HVcbls[length(circ.owp_HVcbls)].elec.volt == 132)
                push!(p132,op.mva)
                push!(w132,op.wnd)
            elseif (circ.owp_HVcbls[length(circ.owp_HVcbls)].elec.volt == 220)
                push!(p220,op.mva)
                push!(w220,op.wnd)
            elseif (circ.owp_HVcbls[length(circ.owp_HVcbls)].elec.volt == 400)
                push!(p400,op.mva)
                push!(w400,op.wnd)
            end
        end
    end
#circ.osss_mog=[]
    #mog transformer
    ossMOG=bus()
    ossMOG.mva=circ.oss_mva
    ossMOG.wnd=circ.oss_wind
    ossMOG.node=opNode
    ossMOG.base_cost=ocn.finance.FC_bld
    ossMOG.num=length(circ.osss_mog)+1
    ossMOG.kV=circ.pcc_cbls[length(circ.pcc_cbls)].elec.volt
    push!(circ.osss_mog,ossMOG)
    circ.osss_mog[length(circ.osss_mog)]=opt_mogXfmrs(circ.osss_mog[length(circ.osss_mog)],pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400,ocn.finance,circ.pcc_cbls[length(circ.pcc_cbls)].elec.volt)

    circ.binary=bn
    circ.decimal=top_bin2dec(bn)
    circ.owpps=buses
    opt_ttlMvCirc(circ)
    return circ
end

#Sums power and wind profiles at connection point
function opt_findWindPwr(buses)
    power_sum=0
    wind_sum=wind()
    wind_sum.pu=zeros(Float32,length(buses[1].wnd.pu))
    wind_sum.ce=zeros(Float32,length(buses[1].wnd.ce))
    wind_sum.delta=0
    wind_sum.lf=0
    for opp in buses
        wind_sum.pu=(wind_sum.pu.+opp.wnd.pu)
        wind_sum.ce=(wind_sum.ce.+opp.wnd.ce)
        wind_sum.delta=(wind_sum.delta+opp.wnd.delta)
        wind_sum.lf=(wind_sum.lf+opp.wnd.lf)
        power_sum=(power_sum+opp.mva)
    end
    wind_sum.pu=(wind_sum.pu)./length(buses)
    wind_sum.ce=(wind_sum.ce)./length(buses)
    wind_sum.delta=(wind_sum.delta)/length(buses)
    wind_sum.lf=(wind_sum.lf)/length(buses)
    return power_sum,wind_sum
end
#**Totals wind profiles of set
function opt_Wsum(wo,wA)
    for w in wA
        wo.pu=(wo.pu.+w.pu)
        wo.ce=(wo.ce.+w.ce)
        wo.delta=(wo.delta+w.delta)
        wo.lf=(wo.lf+w.lf)
    end
    wo.pu=wo.pu./length(wA)
    wo.ce=wo.ce./length(wA)
    wo.delta=wo.delta/length(wA)
    wo.lf=wo.lf/length(wA)
    return wo
end

#**
#places an oss mv_distance towards connection point
function opt_oss4mv2hv(mogNd,owp)
    A=owp.node.xy.y-mogNd.xy.y
    B=mogNd.xy.x-owp.node.xy.x
    C=owp.node.xy.x*mogNd.xy.y-mogNd.xy.x*owp.node.xy.y
    a=A^2+B^2
    if (abs(B)>10^(-3))
        b=2*A*C+2*A*B*owp.node.xy.y-2*B^2*owp.node.xy.x
        c=C^2+2*B*C*owp.node.xy.y-B^2*(owp.zone^2-owp.node.xy.x^2-owp.node.xy.y^2)
        x0=(-b+sqrt(b^2-4*a*c))/(2*a)
        x1=(-b-sqrt(b^2-4*a*c))/(2*a)
        y0=-(A*x0+C)/B
        y1=-(A*x1+C)/B
    else
        b=2*B*C+2*A*B*owp.node.xy.x-2*A^2*owp.node.xy.y
        c=C^2+2*A*C*owp.node.xy.x-A^2*(owp.zone^2-owp.node.xy.x^2-owp.node.xy.y^2)
        y0=(-b+sqrt(b^2-4*a*c))/(2*a)
        y1=(-b-sqrt(b^2-4*a*c))/(2*a)
        x0=-(B*y0+C)/A
        x1=-(B*y1+C)/A
    end
    co_ords=node()
    if ((mogNd.xy.x>=owp.node.xy.x)&&(mogNd.xy.y>owp.node.xy.y))
        #first quadrant
        if ((x0>=owp.node.xy.x)&&(y0>owp.node.xy.y))
            co_ords.xy.x=x0
            co_ords.xy.y=y0
        else
            co_ords.xy.x=x1
            co_ords.xy.y=y1
        end
    elseif ((mogNd.xy.x<owp.node.xy.x)&&(mogNd.xy.y>=owp.node.xy.y))
        #2nd quadrant
        if ((x0<owp.node.xy.x)&&(y0>=owp.node.xy.y))
            co_ords.xy.x=x0
            co_ords.xy.y=y0
        else
            co_ords.xy.x=x1
            co_ords.xy.y=y1
        end
    elseif ((mogNd.xy.x<=owp.node.xy.x)&&(mogNd.xy.y<owp.node.xy.y))
        #3rd quadrant
        if ((x0<=owp.node.xy.x)&&(y0<owp.node.xy.y))
            co_ords.xy.x=x0
            co_ords.xy.y=y0
        else
            co_ords.xy.x=x1
            co_ords.xy.y=y1
        end
    elseif ((mogNd.xy.x>owp.node.xy.x)&&(mogNd.xy.y<=owp.node.xy.y))
        #4th quadrant
        if ((x0>owp.node.xy.x)&&(y0<=owp.node.xy.y))
            co_ords.xy.x=x0
            co_ords.xy.y=y0
        else
            co_ords.xy.x=x1
            co_ords.xy.y=y1
        end
    else
        println("Error: did not match a quandrant in OSS placement at zone radius!!!!")
    end
    return co_ords
end

#binary to decimal conversion
function top_bin2dec(bn)
    dec=0.0
    for (i,bt) in enumerate(bn)
        dec=dec+bt*(2^(i-1))
    end
    return dec
end
############################################################################
#** - not used but may be in the future
function top_dec2bin(dec)
    bn=Int8[]
    while (dec != 0)
        push!(bn,mod(dec,2))
        dec=floor(Int32, dec/2)
    end
    return bn
end
############################################################################

#**
#totals the cost of an HV or MV circuit
function opt_ttlMvCirc(circ)
    circ.cost=0
    for cb in circ.pcc_cbls
        circ.cost=circ.cost+cb.costs.ttl
    end

    for cb in circ.oss2oss_cbls
        circ.cost=circ.cost+cb.costs.ttl
    end

    for cb in circ.owp_MVcbls
        circ.cost=circ.cost+cb.costs.ttl
    end

    for cb in circ.owp_HVcbls
        circ.cost=circ.cost+cb.costs.ttl
    end

    for oss in circ.osss_owp
        circ.cost=circ.cost+oss.base_cost
        for xf in oss.xfmrs
            circ.cost=circ.cost+xf.costs.ttl
        end
    end

    for oss in circ.osss_mog
        circ.cost=circ.cost+oss.base_cost
        for xf in oss.xfmrs
            circ.cost=circ.cost+xf.costs.ttl
        end
    end

    for xf in circ.pcc.xfmrs
        circ.cost=circ.cost+xf.costs.ttl
    end

end
#**
function opt_InitPW()
    pSum_mv=Float32[]
    pSum_hv=Float32[]
    pSum_hv132=Float32[]
    pSum_hv220=Float32[]
    pSum_hv400=Float32[]

    wind_sumMv=wind[]
    wind_sumHv=wind[]
    wind_sum132=wind[]
    wind_sum220=wind[]
    wind_sum400=wind[]
    return pSum_mv,pSum_hv,pSum_hv132,pSum_hv220,pSum_hv400,wind_sumMv,wind_sumHv,wind_sum132,wind_sum220,wind_sum400
end

function opt_combineCblPW(oss_system,circ,pHv,wHv,p132,w132,p220,w220,p400,w400)
    #store wind for mog transformer calcs
    if (oss_system.oss2oss_cbls[length(oss_system.oss2oss_cbls)].elec.volt==oss_system.pcc_cbls[length(oss_system.pcc_cbls)].elec.volt)
        push!(pHv,circ.oss_mva)
        push!(wHv,circ.oss_wind)
    elseif (oss_system.oss2oss_cbls[length(oss_system.oss2oss_cbls)].elec.volt==132)
        push!(p132,circ.oss_mva)
        push!(w132,circ.oss_wind)
    elseif (oss_system.oss2oss_cbls[length(oss_system.oss2oss_cbls)].elec.volt==220)
        push!(p220,circ.oss_mva)
        push!(w220,circ.oss_wind)
    elseif (oss_system.oss2oss_cbls[length(oss_system.oss2oss_cbls)].elec.volt==400)
        push!(p400,circ.oss_mva)
        push!(w400,circ.oss_wind)
    else
    end
    return pHv,wHv,p132,w132,p220,w220,p400,w400
end

#parnt,single_parent
function opt_str8Oss2Oss(circ,oss_system,ocn,pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400)
    #Calculate length and type of cables
    l=lof_pnt2pnt_dist(circ.base_owp.node.xy,oss_system.osss_mog[1].node.xy)
    cb,xs=cstF_MvHvCbloss(l,circ.base_owp.mva,circ.base_owp.wnd,ocn.finance,oss_system.pcc_cbls[1].elec.volt,ocn.sys,ocn.owpps[1].mv_zone,ocn.eqp_data)

    #store cable data
    if (length(cb)==1)
        push!(cb[1].pth,deepcopy(circ.base_owp.node))
        push!(cb[1].pth,deepcopy(oss_system.osss_mog[1].node))
        push!(oss_system.owp_MVcbls,cb[1])
        #store winds and powers
        push!(pMv,circ.base_owp.mva)
        push!(wMv,circ.base_owp.wnd)
    elseif (length(cb)==2)
        #Cable
        mvhvOSSnd=opt_oss4mv2hv(oss_system.osss_mog[1].node,circ.base_owp)
        ocn.buses=ocn.buses+1
        mvhvOSSnd.num=ocn.buses

        push!(cb[1].pth,deepcopy(circ.base_owp.node))
        push!(cb[1].pth,deepcopy(mvhvOSSnd))
        push!(cb[2].pth,deepcopy(mvhvOSSnd))
        push!(cb[2].pth,deepcopy(oss_system.osss_mog[1].node))
        push!(oss_system.owp_MVcbls,deepcopy(cb[1]))
        push!(oss_system.owp_HVcbls,deepcopy(cb[2]))

        #OSS
        ossmv=bus()
        ossmv.node=mvhvOSSnd
        ossmv.base_cost=ocn.finance.FC_bld
        ossmv.mva=circ.base_owp.mva
        ossmv.wnd=circ.base_owp.wnd
        push!(oss_system.osss_owp,ossmv)
        push!(oss_system.osss_owp[length(oss_system.osss_owp)].xfmrs,deepcopy(xs[1]))

        #store winds and powers
        if (oss_system.owp_HVcbls[length(oss_system.owp_HVcbls)].elec.volt == oss_system.pcc_cbls[length(oss_system.pcc_cbls)].elec.volt)
            push!(pHv,circ.base_owp.mva)
            push!(wHv,circ.base_owp.wnd)
        elseif (oss_system.owp_HVcbls[length(oss_system.owp_HVcbls)].elec.volt == 132)
            push!(p132,circ.base_owp.mva)
            push!(w132,circ.base_owp.wnd)
        elseif (oss_system.owp_HVcbls[length(oss_system.owp_HVcbls)].elec.volt == 220)
            push!(p220,circ.base_owp.mva)
            push!(w220,circ.base_owp.wnd)
        elseif (oss_system.owp_HVcbls[length(oss_system.owp_HVcbls)].elec.volt == 400)
            push!(p400,circ.base_owp.mva)
            push!(w400,circ.base_owp.wnd)
        end
    end
    return oss_system,pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400
end
