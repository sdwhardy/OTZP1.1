
#=
ocn=ocean
pcc=ocn.pccs[2]
owpps=ocn.owpps
paths=opt_mvOSSplacement(ocn,owpps,pcc)
=#
function opt_compoundOSS(ocn,pcc)
    mhv_circs=ocn.circuits
    owpps_tbl_mhv=Array{Array{Int8, 1}, 1}()
    oss_system=circuit()
    oss_QFl=Array{circuit,1}()
    oss_QWrk=Array{circuit,1}()
    rez=Array{circuit,1}()
    term_Node=ocn.owpps[1].num
    for crc in mhv_circs
        if ((crc.owpps[1].num!=term_Node) && (sum(crc.binary)!=1))
            crc.Qing=true
            push!(oss_QFl,crc)
        end
    end
    indx2=1
    oss_QWrk=oss_QFl[indx2]
    while (oss_QWrk.Qing != false)
        parent_Osss=findall(x->x==1,oss_QWrk.binary)
        least_parent=parent_Osss[1]
        child_Osss=findall(x->x==0, oss_QWrk.binary)
        homeless_child=findall(x->x>least_parent, child_Osss)
        homeowning_child=findall(x->x<least_parent, child_Osss)
        clms=trunc(Int,length(child_Osss)+1)
        rows=trunc(Int, 2.0^clms)
        empty_tbl_dumb=eensF_blankTbl(rows,clms)
        empty_tbl=Array{Array{Int8, 1}, 1}()
        for indx=1:length(empty_tbl_dumb[:,1])
            if ((sum(empty_tbl_dumb[indx,1:length(homeowning_child)])>0)&&(empty_tbl_dumb[indx,least_parent]==1))
                push!(empty_tbl,empty_tbl_dumb[indx,:])
            end
        end
        for indx0=1:length(empty_tbl[:,1])
        #indx0=1
            bs2add=bus[]
            bn=deepcopy(oss_QWrk.binary)
            for hoc in homeowning_child
                if (empty_tbl[indx0][hoc]==1)
                    bn[hoc]=1
                    push!(bs2add,ocn.owpps[hoc])
                end
            end
            for hlc in homeless_child
                if (empty_tbl[indx0][hlc+1]==1)
                    bn[child_Osss[hlc]]=1
                    push!(bs2add,ocn.owpps[child_Osss[hlc]])
                end
            end

            oss_system, domain_oss,domain_edges=opt_compoundOssSystem(bs2add,pcc,oss_QWrk,bn,ocn)
            oss_QWrk.Qing=false
            if (mhv_circs[oss_system.decimal].cost > oss_system.cost)
                mhv_circs[oss_system.decimal]=deepcopy(oss_system)
                ocn.discretedom.nodes=deepcopy(domain_oss)
                ocn.discretedom.edges=deepcopy(domain_edges)
                if (oss_system.binary[1]==0)
                    indx1=1
                    while (oss_system.decimal != oss_QFl[indx1].decimal)
                        indx1=indx1+1
                    end
                    oss_system.Qing=true
                    oss_QFl[indx1]=deepcopy(oss_system)
                end
            else
            end
            #=push!(rez,deepcopy(oss_system))
            ocn.discretedom.nodes=deepcopy(domain_oss)
            ocn.discretedom.edges=deepcopy(domain_edges)
            if (oss_system.binary[1]==0)
                indx1=1
                while (oss_system.decimal != oss_QFl[indx1].decimal)
                    indx1=indx1+1
                end
                oss_system.Qing=true
                #println("oss_QFl[indx1] b4-1: "*string(length(oss_QFl[indx1].oss_wind.pu)))
                oss_QFl[indx1]=deepcopy(oss_system)
                #println("oss_QFl[indx1] b4-2: "*string(length(oss_QFl[indx1].oss_wind.pu)))
            end=#
        end
        oss_QWrk.Qing = false
        while (indx2 <= length(oss_QFl) && oss_QFl[indx2].Qing==false)
            indx2=indx2+1
        end
        if (indx2 <= length(oss_QFl))
            oss_QWrk=oss_QFl[indx2]
            println("position in oss_QFl: "*string(indx2)*"/"*string(length(oss_QFl)))
        end
    end
    return mhv_circs
end
#=
buses=bs2add
parnt=oss_QWrk
=#
function opt_compoundOssSystem(buses,pcc,parnt,bn,ocn)
    circ=circuit()
    circ.binary=bn
    circ.parent_circ=parnt.decimal
    circ.decimal=top_bin2dec(bn)
    circ.pcc=deepcopy(pcc)
    for (indx,owp) in enumerate(bn)
        if (owp==1)
            push!(circ.owpps,ocn.owpps[indx])
        end
    end
    owp=buses[1]
    oss_node,domain_oss,edges_oss,power_sum,wind_sum=opt_mvhvOss1stLocal(owp,buses,pcc,ocn)

    #xys=opt_mvhvOss1stLocal(owp,buses,pcc,ocn)
    oss=bus()
    oss.node=deepcopy(parnt.osss_mog[length(parnt.osss_mog)].node)
    oss.wnd=deepcopy(parnt.osss_mog[length(parnt.osss_mog)].wnd)
    oss.mva=deepcopy(parnt.osss_mog[length(parnt.osss_mog)].mva)
    #push!(xys,oss.node.xy)
    push!(buses,oss)

    #domain_oss,oss_node,power_sum,wind_sum,edges_oss=opt_adjustPath(owp,ocn,xys,buses,pcc)
    push!(circ.pths,deepcopy(as_Astar(domain_oss[pcc.node.num],oss_node,domain_oss)))
    push!(circ.lengths,circ.pths[length(circ.pths)].G_cost)
    cbl_xfo=cstF_HvCblallKvo2p(circ.lengths[length(circ.lengths)],power_sum,wind_sum,ocn.finance,pcc)
    #here
    currentNode=circ.pths[1]
    goalNode=circ.pths[1].goal
    push!(cbl_xfo[1].pth,deepcopy(currentNode))
    while currentNode.num != goalNode
        push!(cbl_xfo[1].pth,deepcopy(currentNode.parent))
        currentNode=currentNode.parent
    end
    push!(circ.pcc_cbls,deepcopy(cbl_xfo[1]))
    circ.osss_owp=deepcopy(parnt.osss_owp)
    circ.osss_mog=deepcopy(parnt.osss_mog)
    circ.owp_MVcbls=deepcopy(parnt.owp_MVcbls)
    circ.owp_HVcbls=deepcopy(parnt.owp_HVcbls)
    circ.oss2oss_cbls=deepcopy(parnt.oss2oss_cbls)
    push!(circ.pths,deepcopy(as_Astar(domain_oss[parnt.osss_mog[length(parnt.osss_mog)].node.num],oss_node,domain_oss)))
    push!(circ.oss2oss_cbls,deepcopy(cstF_HvCblo2o(circ.pths[length(circ.pths)].G_cost,parnt.osss_mog[length(parnt.osss_mog)].mva,parnt.pcc_cbls[length(parnt.pcc_cbls)].elec.volt,parnt.osss_mog[length(parnt.osss_mog)].wnd,ocn.finance)))
    currentNode=circ.pths[length(circ.pths)]
    goalNode=circ.pths[length(circ.pths)].goal
    push!(circ.oss2oss_cbls[length(circ.oss2oss_cbls)].pth,deepcopy(currentNode))
    while currentNode.num != goalNode
        push!(circ.oss2oss_cbls[length(circ.oss2oss_cbls)].pth,deepcopy(currentNode.parent))
        currentNode=currentNode.parent
    end
    circ.oss2oss_cbls[length(circ.oss2oss_cbls)].pth=reverse!(circ.oss2oss_cbls[length(circ.oss2oss_cbls)].pth)
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
    if (circ.oss2oss_cbls[length(circ.oss2oss_cbls)].elec.volt==circ.pcc_cbls[length(circ.pcc_cbls)].elec.volt)
        iSum_hv=iSum_hv+1
        pSum_hv=pSum_hv+10
        wind_sumHv.pu=(wind_sumHv.pu.+deepcopy(parnt.osss_mog[length(parnt.osss_mog)].wnd.pu))
        wind_sumHv.ce=(wind_sumHv.ce.+deepcopy(parnt.osss_mog[length(parnt.osss_mog)].wnd.ce))
        wind_sumHv.delta=(wind_sumHv.delta+deepcopy(parnt.osss_mog[length(parnt.osss_mog)].wnd.delta))
        wind_sumHv.lf=(wind_sumHv.lf+deepcopy(parnt.osss_mog[length(parnt.osss_mog)].wnd.lf))
    elseif (circ.oss2oss_cbls[length(circ.oss2oss_cbls)].elec.volt==132)
        iSum_hv132=iSum_hv132+1
        pSum_hv132=pSum_hv132+deepcopy(parnt.osss_mog[length(parnt.osss_mog)].mva)
        wind_sum132.pu=(wind_sum132.pu.+deepcopy(parnt.osss_mog[length(parnt.osss_mog)].wnd.pu))
        wind_sum132.ce=(wind_sum132.ce.+deepcopy(parnt.osss_mog[length(parnt.osss_mog)].wnd.ce))
        wind_sum132.delta=(wind_sum132.delta+deepcopy(parnt.osss_mog[length(parnt.osss_mog)].wnd.delta))
        wind_sum132.lf=(wind_sum132.lf+deepcopy(parnt.osss_mog[length(parnt.osss_mog)].wnd.lf))
    elseif (circ.oss2oss_cbls[length(circ.oss2oss_cbls)].elec.volt==220)
        iSum_hv220=iSum_hv220+1
        pSum_hv220=pSum_hv220+deepcopy(parnt.osss_mog[length(parnt.osss_mog)].mva)
        wind_sum220.pu=(wind_sum220.pu.+deepcopy(parnt.osss_mog[length(parnt.osss_mog)].wnd.pu))
        wind_sum220.ce=(wind_sum220.ce.+deepcopy(parnt.osss_mog[length(parnt.osss_mog)].wnd.ce))
        wind_sum220.delta=(wind_sum220.delta+deepcopy(parnt.osss_mog[length(parnt.osss_mog)].wnd.delta))
        wind_sum220.lf=(wind_sum220.lf+deepcopy(parnt.osss_mog[length(parnt.osss_mog)].wnd.lf))
    elseif (circ.oss2oss_cbls[length(circ.oss2oss_cbls)].elec.volt==400)
        iSum_hv400=iSum_hv400+1
        pSum_hv400=pSum_hv400+deepcopy(parnt.osss_mog[length(parnt.osss_mog)].mva)
        wind_sum400.pu=(wind_sum400.pu.+deepcopy(parnt.osss_mog[length(parnt.osss_mog)].wnd.pu))
        wind_sum400.ce=(wind_sum400.ce.+deepcopy(parnt.osss_mog[length(parnt.osss_mog)].wnd.ce))
        wind_sum400.delta=(wind_sum400.delta+deepcopy(parnt.osss_mog[length(parnt.osss_mog)].wnd.delta))
        wind_sum400.lf=(wind_sum400.lf+deepcopy(parnt.osss_mog[length(parnt.osss_mog)].wnd.lf))
    else
    end
    for owp in buses[1:length(buses)-1]
        push!(circ.pths,deepcopy(as_Astar(domain_oss[owp.node.num],oss_node,domain_oss)))
        cs,xs=cstF_MvHvCbloss(circ.pths[length(circ.pths)].G_cost,owp.mva,owp.wnd,ocn.finance,pcc,cbl_xfo[1].elec.volt)
        currentNode=circ.pths[length(circ.pths)]
        goalNode=circ.pths[length(circ.pths)].goal
        if length(cs) == 1
            while currentNode.num != goalNode
                push!(cs[1].pth,deepcopy(currentNode))
                currentNode=currentNode.parent
            end
            push!(cs[1].pth,deepcopy(currentNode))
            cs[1].pth=reverse!(cs[1].pth)
            push!(circ.owp_MVcbls,deepcopy(cs[1]))
        elseif length(cs) == 2
            push!(cs[2].pth,deepcopy(currentNode))
            while currentNode.parent.num != goalNode
                push!(cs[2].pth,deepcopy(currentNode.parent))
                currentNode=currentNode.parent
            end
            push!(cs[1].pth,deepcopy(currentNode))
            push!(cs[1].pth,deepcopy(currentNode.parent))
            cs[1].pth=reverse!(cs[1].pth)
            cs[2].pth=reverse!(cs[2].pth)
            push!(circ.owp_MVcbls,deepcopy(cs[1]))
            push!(circ.owp_HVcbls,deepcopy(cs[2]))
        end
        if length(cs) == 1
            iSum_mv=iSum_mv+1
            pSum_mv=pSum_mv+deepcopy(owp.mva)
            wind_sumMv.pu=deepcopy((wind_sumMv.pu.+owp.wnd.pu))
            wind_sumMv.ce=deepcopy((wind_sumMv.ce.+owp.wnd.ce))
            wind_sumMv.delta=deepcopy((wind_sumMv.delta+owp.wnd.delta))
            wind_sumMv.lf=deepcopy((wind_sumMv.lf+owp.wnd.lf))
        elseif cs[2].elec.volt == cbl_xfo[1].elec.volt
            ossBus=bus()
            ossBus.mva=owp.mva
            ossBus.wnd=owp.wnd
            ossBus.node=circ.owp_MVcbls[length(circ.owp_MVcbls)].pth[length(circ.owp_MVcbls[length(circ.owp_MVcbls)].pth)]
            push!(ossBus.xfmrs,deepcopy(xs[1]))
            push!(circ.osss_owp,deepcopy(ossBus))
            iSum_hv=iSum_hv+1
            pSum_hv=pSum_hv+10
            wind_sumHv.pu=(wind_sumHv.pu.+deepcopy(owp.wnd.pu))
            wind_sumHv.ce=(wind_sumHv.ce.+deepcopy(owp.wnd.ce))
            wind_sumHv.delta=(wind_sumHv.delta+deepcopy(owp.wnd.delta))
            wind_sumHv.lf=(wind_sumHv.lf+deepcopy(owp.wnd.lf))
        elseif cs[2].elec.volt == 132
            ossBus=bus()
            ossBus.mva=owp.mva
            ossBus.wnd=owp.wnd
            ossBus.node=circ.owp_MVcbls[length(circ.owp_MVcbls)].pth[length(circ.owp_MVcbls[length(circ.owp_MVcbls)].pth)]
            push!(ossBus.xfmrs,deepcopy(xs[1]))
            push!(circ.osss_owp,deepcopy(ossBus))
            iSum_hv132=iSum_hv132+1
            pSum_hv132=pSum_hv132+deepcopy(owp.mva)
            wind_sum132.pu=(wind_sum132.pu.+deepcopy(owp.wnd.pu))
            wind_sum132.ce=(wind_sum132.ce.+deepcopy(owp.wnd.ce))
            wind_sum132.delta=(wind_sum132.delta+deepcopy(owp.wnd.delta))
            wind_sum132.lf=(wind_sum132.lf+deepcopy(owp.wnd.lf))
        elseif cs[2].elec.volt == 220
            ossBus=bus()
            ossBus.mva=owp.mva
            ossBus.wnd=owp.wnd
            ossBus.node=circ.owp_MVcbls[length(circ.owp_MVcbls)].pth[length(circ.owp_MVcbls[length(circ.owp_MVcbls)].pth)]
            push!(ossBus.xfmrs,deepcopy(xs[1]))
            push!(circ.osss_owp,deepcopy(ossBus))
            iSum_hv220=iSum_hv220+1
            pSum_hv220=pSum_hv220+deepcopy(owp.mva)
            wind_sum220.pu=(wind_sum220.pu.+deepcopy(owp.wnd.pu))
            wind_sum220.ce=(wind_sum220.ce.+deepcopy(owp.wnd.ce))
            wind_sum220.delta=(wind_sum220.delta+deepcopy(owp.wnd.delta))
            wind_sum220.lf=(wind_sum220.lf+deepcopy(owp.wnd.lf))
        elseif cs[2].elec.volt == 400
            ossBus=bus()
            ossBus.mva=owp.mva
            ossBus.wnd=owp.wnd
            ossBus.node=circ.owp_MVcbls[length(circ.owp_MVcbls)].pth[length(circ.owp_MVcbls[length(circ.owp_MVcbls)].pth)]
            push!(ossBus.xfmrs,deepcopy(xs[1]))
            push!(circ.osss_owp,deepcopy(ossBus))
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
        ossMOG=bus()
        pSum_hv=pSum_hv+100
        ossMOG.node=oss_node
        push!(ossMOG.xfmrs,cstF_xfo_oss(pSum_hv,wind_sumHv,ocn.finance))
        push!(circ.osss_mog,deepcopy(ossMOG))
    elseif pSum_hv != 0
        ossMOG=bus()
        ossMOG.node=oss_node
        push!(ossMOG.xfmrs,cstF_xfo_oss(pSum_hv,wind_sumHv,ocn.finance))
        push!(circ.osss_mog,deepcopy(ossMOG))
    else
        ossMOG=bus()
        ossMOG.node=oss_node
        push!(circ.osss_mog,ossMOG)
    end
    if pSum_mv != 0
        push!(circ.osss_mog[length(circ.osss_mog)].xfmrs,deepcopy(cstF_xfo_oss(pSum_mv,wind_sumMv,ocn.finance)))
    end
    if pSum_hv132 != 0
        push!(circ.osss_mog[length(circ.osss_mog)].xfmrs,deepcopy(cstF_xfo_oss(pSum_hv132,wind_sum132,ocn.finance)))
    end
    if pSum_hv220 != 0
        push!(circ.osss_mog[length(circ.osss_mog)].xfmrs,deepcopy(cstF_xfo_oss(pSum_hv220,wind_sum220,ocn.finance)))
    end
    if pSum_hv400 != 0
        push!(circ.osss_mog[length(circ.osss_mog)].xfmrs,deepcopy(cstF_xfo_oss(pSum_hv400,wind_sum400,ocn.finance)))
    end
    circ.osss_mog[length(circ.osss_mog)].mva=power_sum
    circ.osss_mog[length(circ.osss_mog)].wnd=deepcopy(wind_sum)
    circ.pcc=deepcopy(pcc)
    if (pcc.kV != circ.pcc_cbls[1].elec.volt)
        println("power: "*string(power_sum))
        push!(circ.pcc.xfmrs,deepcopy(cstF_xfo_pcc(power_sum,wind_sum,ocn.finance)))
    end

    circ.binary=bn
    circ.decimal=top_bin2dec(bn)
    circ.base_owp=owp
    opt_ttlMvCirc(circ)
    #println("circ return: "*string(length(circ.osss_mog[length(circ.osss_mog)].wnd.pu)))
    return circ, domain_oss,edges_oss
end


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
    #for indx0=20:20
        bus_dummies=bus[]
        for (indx1,owp) in enumerate(owpps_tbl_hv[indx0])
            if (owp==1)
                push!(bus_dummies,ocn.owpps[indx1])
            end
        end
        oss_system, discrete_dom, discrete_edges=opt_hvOssSystem(bus_dummies,pcc,ocn,owpps_tbl_hv[indx0])
        if (oss_system.cost<mv_circs[indx0].cost)
            mv_circs[indx0]=deepcopy(oss_system)
            ocn.discretedom.nodes=deepcopy(discrete_dom)
            ocn.discretedom.edges=deepcopy(discrete_edges)
        end
    end
    return mv_circs
end

function opt_keepBestMvHv(mv_cs,hv_cs)
    oss_stms=Array{circuit,1}()
    for mv_c in mv_cs
        for hv_c in hv_cs
            if (mv_c.decimal==hv_c.decimal)
                if (mv_c.cost <=  hv_c.cost)
                    push!(oss_stms,deepcopy(mv_c))
                else
                    push!(oss_stms,deepcopy(hv_c))
                end
            else
            end
        end
    end
    return oss_stms
end
#length(owpps_tbl[:,1])
#length(owpps_tbl_hv[:,1])
function opt_mvOSSplacement(ocn,owpps,pcc)
    #owpps is an array of owpps ordered closest to farthest from the designated pcc
    owpps_tbl=top_mvTopos(ocn.owpps)
    owpps_tbl=owpps_tbl[end:-1:1,end:-1:1]
    bus_dummies=Array{bus,1}()
    oss_system=circuit()
    oss_systems=Array{circuit,1}()
    for indx0=1:length(owpps_tbl[:,1])
    #for indx0=4:4
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
    circ.pcc=deepcopy(pcc)
    push!(circ.owpps,owp)
    push!(circ.pths,as_Astar(owp.node,pcc.node,ocn.discretedom.nodes))
    push!(circ.lengths,circ.pths[1].G_cost)
    #not necessarily M
    cs=cstF_MvHvCblpcc(circ.pths[1].G_cost,owp.mva,owp.wnd,ocn.finance,pcc)
    currentNode=circ.pths[1]
    goalNode=circ.pths[1].goal
    if length(cs) == 1
        while currentNode.parent.num != goalNode
            push!(cs[1].pth,deepcopy(currentNode))
            currentNode=currentNode.parent
        end
        push!(circ.owp_MVcbls,deepcopy(cs[1]))
    elseif length(cs) == 2
        push!(cs[2].pth,deepcopy(currentNode))
        while currentNode.parent.num != goalNode
            push!(cs[2].pth,deepcopy(currentNode.parent))
            currentNode=currentNode.parent
        end
        push!(cs[1].pth,deepcopy(currentNode))
        push!(cs[1].pth,deepcopy(currentNode.parent))
        push!(circ.owp_MVcbls,deepcopy(cs[1]))
        push!(circ.pcc_cbls,deepcopy(cs[2]))
    end
    circ.base_owp=owp

    if (length(cs) > 1)
        ossmv=bus()
        ossmv.node=currentNode
        push!(circ.osss_owp,ossmv)
        push!(circ.osss_owp[1].xfmrs,deepcopy(cstF_xfo_oss(owp.mva,owp.wnd,ocn.finance)))
    end
    pccBus=bus()
    pccBus.node=circ.pths[1]
    circ.pcc=pccBus
    if (cs[length(cs)].elec.volt != pcc.kV)
        push!(circ.pcc.xfmrs,deepcopy(cstF_xfo_pcc(owp.mva,owp.wnd,ocn.finance)))
    end
    opt_ttlMvCirc(circ)

    return circ
end
#buses=bus_dummies
#bn=owpps_tbl_hv[indx0]

function opt_hvOssSystem(buses,pcc,ocn,bn)
    owp=buses[1]
    #domain_oss,oss_node,power_sum,wind_sum,edges_oss=opt_adjustPath(owp,ocn,xys,buses,pcc)
    oss_node,domain_oss,edges_oss,power_sum,wind_sum=opt_mvhvOss1stLocal(owp,buses,pcc,ocn)
    circ=circuit()
    circ.oss_wind=deepcopy(wind_sum)
    circ.oss_mva=deepcopy(power_sum)
    push!(circ.pths,deepcopy(as_Astar(domain_oss[pcc.node.num],oss_node,domain_oss)))
    push!(circ.lengths,circ.pths[length(circ.pths)].G_cost)
    println("mva: "*string(power_sum))
    cbl_xfo=cstF_HvCblallKvo2p(circ.lengths[length(circ.lengths)],power_sum,wind_sum,ocn.finance,pcc)
    currentNode=circ.pths[length(circ.pths)]
    goalNode=circ.pths[length(circ.pths)].goal
    push!(cbl_xfo[1].pth,deepcopy(currentNode))
    while currentNode.num != goalNode
        push!(cbl_xfo[1].pth,deepcopy(currentNode.parent))
        currentNode=currentNode.parent
    end
    push!(circ.pcc_cbls,cbl_xfo[1])
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
        cs,xs=cstF_MvHvCbloss(circ.pths[length(circ.pths)].G_cost,owp.mva,owp.wnd,ocn.finance,pcc,cbl_xfo[1].elec.volt)
        currentNode=circ.pths[length(circ.pths)]
        goalNode=circ.pths[length(circ.pths)].goal
        if length(cs) == 1
            while currentNode.num != goalNode
                push!(cs[1].pth,deepcopy(currentNode))
                currentNode=currentNode.parent
            end
            push!(cs[1].pth,deepcopy(currentNode))
            push!(circ.owp_MVcbls,deepcopy(cs[1]))
        elseif length(cs) == 2
            push!(cs[2].pth,deepcopy(currentNode))
            while currentNode.parent.num != goalNode
                push!(cs[2].pth,deepcopy(currentNode.parent))
                currentNode=currentNode.parent
            end
            push!(cs[1].pth,deepcopy(currentNode))
            push!(cs[1].pth,deepcopy(currentNode.parent))
            push!(circ.owp_MVcbls,deepcopy(cs[1]))
            push!(circ.owp_HVcbls,deepcopy(cs[2]))
        end
        if length(cs) == 1
            iSum_mv=iSum_mv+1
            pSum_mv=pSum_mv+deepcopy(owp.mva)
            wind_sumMv.pu=deepcopy((wind_sumMv.pu.+owp.wnd.pu))
            wind_sumMv.ce=deepcopy((wind_sumMv.ce.+owp.wnd.ce))
            wind_sumMv.delta=deepcopy((wind_sumMv.delta+owp.wnd.delta))
            wind_sumMv.lf=deepcopy((wind_sumMv.lf+owp.wnd.lf))
        elseif cs[2].elec.volt == cbl_xfo[1].elec.volt
            ossBus=bus()
            ossBus.mva=owp.mva
            ossBus.wnd=owp.wnd
            ossBus.node=cs[2].pth[length(cs[2].pth)]
            push!(ossBus.xfmrs,xs[1])
            push!(circ.osss_owp,ossBus)
            iSum_hv=iSum_hv+1
            pSum_hv=pSum_hv+10
            wind_sumHv.pu=(wind_sumHv.pu.+deepcopy(owp.wnd.pu))
            wind_sumHv.ce=(wind_sumHv.ce.+deepcopy(owp.wnd.ce))
            wind_sumHv.delta=(wind_sumHv.delta+deepcopy(owp.wnd.delta))
            wind_sumHv.lf=(wind_sumHv.lf+deepcopy(owp.wnd.lf))
        elseif cs[2].elec.volt == 132
            ossBus=bus()
            ossBus.mva=owp.mva
            ossBus.wnd=owp.wnd
            ossBus.node=cs[2].pth[length(cs[2].pth)]
            push!(ossBus.xfmrs,xs[1])
            push!(circ.osss_owp,ossBus)
            iSum_hv132=iSum_hv132+1
            pSum_hv132=pSum_hv132+deepcopy(owp.mva)
            wind_sum132.pu=(wind_sum132.pu.+deepcopy(owp.wnd.pu))
            wind_sum132.ce=(wind_sum132.ce.+deepcopy(owp.wnd.ce))
            wind_sum132.delta=(wind_sum132.delta+deepcopy(owp.wnd.delta))
            wind_sum132.lf=(wind_sum132.lf+deepcopy(owp.wnd.lf))
        elseif cs[2].elec.volt == 220
            ossBus=bus()
            ossBus.mva=owp.mva
            ossBus.wnd=owp.wnd
            ossBus.node=cs[2].pth[length(cs[2].pth)]
            push!(ossBus.xfmrs,xs[1])
            push!(circ.osss_owp,ossBus)
            iSum_hv220=iSum_hv220+1
            pSum_hv220=pSum_hv220+deepcopy(owp.mva)
            wind_sum220.pu=(wind_sum220.pu.+deepcopy(owp.wnd.pu))
            wind_sum220.ce=(wind_sum220.ce.+deepcopy(owp.wnd.ce))
            wind_sum220.delta=(wind_sum220.delta+deepcopy(owp.wnd.delta))
            wind_sum220.lf=(wind_sum220.lf+deepcopy(owp.wnd.lf))
        elseif cs[2].elec.volt == 400
            ossBus=bus()
            ossBus.mva=owp.mva
            ossBus.wnd=owp.wnd
            ossBus.node=cs[2].pth[length(cs[2].pth)]
            push!(ossBus.xfmrs,xs[1])
            push!(circ.osss_owp,ossBus)
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
        ossMOG=bus()
        pSum_hv=pSum_hv+100
        ossMOG.node=oss_node
        push!(circ.osss_mog,ossMOG)
        push!(circ.osss_mog[length(circ.osss_mog)].xfmrs,cstF_xfo_oss(pSum_hv,wind_sumHv,ocn.finance))
    elseif pSum_hv != 0
        ossMOG=bus()
        ossMOG.node=oss_node
        push!(circ.osss_mog,ossMOG)
        push!(circ.osss_mog[length(circ.osss_mog)].xfmrs,cstF_xfo_oss(pSum_hv,wind_sumHv,ocn.finance))
    else
        ossMOG=bus()
        ossMOG.node=oss_node
        push!(circ.osss_mog,ossMOG)
    end
    if pSum_mv != 0
        push!(circ.osss_mog[length(circ.osss_mog)].xfmrs,cstF_xfo_oss(pSum_mv,wind_sumMv,ocn.finance))
    end
    if pSum_hv132 != 0
        push!(circ.osss_mog[length(circ.osss_mog)].xfmrs,cstF_xfo_oss(pSum_hv132,wind_sum132,ocn.finance))
    end
    if pSum_hv220 != 0
        push!(circ.osss_mog[length(circ.osss_mog)].xfmrs,cstF_xfo_oss(pSum_hv220,wind_sum220,ocn.finance))
    end
    if pSum_hv400 != 0
        push!(circ.osss_mog[length(circ.osss_mog)].xfmrs,cstF_xfo_oss(pSum_hv400,wind_sum400,ocn.finance))
    end
    ossMOG.mva=power_sum
    ossMOG.wnd=wind_sum
    circ.pcc=deepcopy(pcc)

    if (pcc.kV != circ.pcc_cbls[1].elec.volt)
        println("power: "*string(power_sum))
        push!(circ.pcc.xfmrs,cstF_xfo_pcc(power_sum,wind_sum,ocn.finance))
    end

    circ.binary=bn
    circ.decimal=top_bin2dec(bn)
    circ.owpps=buses
    circ.base_owp=owp
    opt_ttlMvCirc(circ)
    return circ,domain_oss,edges_oss
end

#bn=owpps_tbl[indx0,:]
#buses=bus_dummies
function opt_mvOssSystem(mv_square,buses,pcc,ocn,bn)
    owp=opt_buses2nodes4MVoptLocal(mv_square,ocn.owpps,buses)
    oss_node,domain_oss,edges_oss,power_sum,wind_sum=opt_mvhvOss1stLocal(owp,buses,pcc,ocn)
    #domain_oss,oss_node,power_sum,wind_sum,edges_oss=opt_adjustPath(owp,ocn,xys,buses,pcc)
    circ=circuit()
    circ.oss_wind=deepcopy(wind_sum)
    circ.oss_mva=deepcopy(power_sum)
    push!(circ.pths,deepcopy(as_Astar(domain_oss[pcc.node.num],oss_node,domain_oss)))
    push!(circ.lengths,circ.pths[length(circ.pths)].G_cost)
    #################### This needs to be a high voltage cable!!!
    println("mva: "*string(power_sum))
    cbl_xfo=cstF_HvCblallKvo2p(circ.lengths[length(circ.lengths)],power_sum,wind_sum,ocn.finance,pcc)
    currentNode=circ.pths[1]
    goalNode=circ.pths[1].goal
    push!(cbl_xfo[1].pth,deepcopy(currentNode))
    while currentNode.num != goalNode
        push!(cbl_xfo[1].pth,deepcopy(currentNode.parent))
        currentNode=currentNode.parent
    end
    push!(circ.pcc_cbls,cbl_xfo[1])
    ossMOG=bus()
    ossMOG.mva=power_sum
    ossMOG.wnd=wind_sum
    ossMOG.node=oss_node
    push!(circ.osss_mog,ossMOG)
    push!(circ.osss_mog[1].xfmrs,cstF_xfo_oss(power_sum,wind_sum,ocn.finance))
    circ.pcc=deepcopy(pcc)
    if (circ.pcc_cbls[length(circ.pcc_cbls)].elec.volt != circ.pcc.kV)
        push!(circ.pcc.xfmrs,cbl_xfo[2])
    end
    for owp in buses
        push!(circ.pths,deepcopy(as_Astar(domain_oss[owp.node.num],oss_node,domain_oss)))
        push!(circ.lengths,circ.pths[length(circ.pths)].G_cost)
        cb=cstF_MvCbl3366(circ.lengths[length(circ.lengths)],owp.mva,owp.wnd,ocn.finance)
        currentNode=circ.pths[length(circ.pths)]
        goalNode=circ.pths[length(circ.pths)].goal
        push!(cb.pth,deepcopy(currentNode))
        while currentNode.num != goalNode
            push!(cb.pth,deepcopy(currentNode.parent))
            currentNode=currentNode.parent
        end
        push!(circ.owp_MVcbls,cb)
    end


    circ.binary=bn
    circ.decimal=top_bin2dec(bn)
    #circ.pcc_length=69.69
    circ.owpps=buses
    #circ.owp_lengths=69.69
    circ.base_owp=owp
    opt_ttlMvCirc(circ)
    ocn.discretedom.nodes=deepcopy(domain_oss)
    ocn.discretedom.edges=deepcopy(edges_oss)
    return circ
end

function opt_MakeConstraints(xys,owp,owps,pcc)
#=    #find major axis
    ymajor=true
    if (abs(pcc.node.xy.y-owp.node.xy.y)>=abs(pcc.node.xy.x-owp.node.xy.x))
        ymajor=true
    else
        ymajor=false
    end

    #find minor axis mean
    minor_ave=0
    for xy in xys[2:length(xys)]
        if ymajor==true
            minor_ave=minor_ave+deepcopy(xy.x)
        else
            minor_ave=minor_ave+deepcopy(xy.y)
        end
    end
    minor_ave=minor_ave/length(xys[2:length(xys)])

    #lower/upper lim
    if (ymajor==true)
        cntr=(2*owp.node.xy.x-pcc.node.xy.x-minor_ave)
        #western connection
        if (cntr>=0)
            upper_lim=owps[1].node.xy.x-owps[1].farm.neg_width
            for opp in owps[1:owp.num]
                if opp.
            end
        end

    if (ymajor==true) && (<=0))
        lwr_lim=#
end

function opt_adjustPath(owp,ocn,xys,buses,pcc)
    paths=node[]
    opt_MakeConstraints(xys,owp,ocn.owpps,pcc)
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
    #error
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
    for (i,owpp) in enumerate(buses)
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
    return domain_oss,oss_node,power_sum,wind_sum,edges_oss
end

function opt_mvhvOss1stLocal(owp,buses,pcc,ocn)
    #find PCC major axis
    ymajorPCC=true
    if (abs(pcc.node.xy.y-owp.node.xy.y)>=abs(pcc.node.xy.x-owp.node.xy.x))
        ymajorPCC=true
    else
        ymajorPCC=false
    end

    ymajorOWPP=true
    if (abs(buses[1].node.xy.y-buses[length(buses)].node.xy.y)>=abs(buses[1].node.xy.x-buses[length(buses)].node.xy.x))
        ymajorOWPP=true
    else
        ymajorOWPP=false
    end

    #find min and max connection points
    lwr_xy=xy()
    upr_xy=xy()
    cnt_xy=xy()
    if (ymajorPCC==true)
        if (pcc.node.xy.y>=owp.node.xy.y)
            lwr_xy.y=owp.node.xy.y+owp.zone.pos_height
        else
            lwr_xy.y=owp.node.xy.y-owp.zone.neg_height
        end
        upr_xy.y=lwr_xy.y
        cnt_xy.y=lwr_xy.y
        lwr_xy.x=owp.node.xy.x-owp.zone.neg_width
        upr_xy.x=owp.node.xy.x+owp.zone.pos_width
        if (owp.node.xy.x-owp.zone.neg_width<pcc.node.xy.x && pcc.node.xy.x<owp.node.xy.x+owp.zone.pos_width)
            cnt_xy.x=pcc.node.xy.x
        else
            cnt_xy.x=owp.node.xy.x
        end
    else
        if (pcc.node.xy.x>=owp.node.xy.y)
            lwr_xy.x=owp.node.xy.x+owp.zone.pos_width
        else
            lwr_xy.x=owp.node.xy.x-owp.zone.neg_width
        end
        upr_xy.x=lwr_xy.x
        cnt_xy.x=lwr_xy.x
        lwr_xy.y=owp.node.xy.y-owp.zone.neg_height
        upr_xy.y=owp.node.xy.y+owp.zone.pos_height

        if (owp.node.xy.y-owp.zone.neg_height<pcc.node.xy.y && pcc.node.xy.y<owp.node.xy.y+owp.zone.pos_height)
            cnt_xy.y=pcc.node.xy.y
        else
            cnt_xy.y=owp.node.xy.y
        end
    end

    #find nodes
    lwr_node=node()
    upr_node=node()
    cnt_node=node()
    lwr_node_mag=Inf
    upr_node_mag=Inf
    cnt_node_mag=Inf

    for nd in (owp.zone.nodes)
        tl_mag=lof_pnt2pnt_dist(nd.xy,lwr_xy)
        if (tl_mag<lwr_node_mag)
            lwr_node_mag=deepcopy(tl_mag)
            lwr_node=deepcopy(nd)
        end
        tu_mag=lof_pnt2pnt_dist(nd.xy,upr_xy)
        if (tu_mag<upr_node_mag)
            upr_node_mag=deepcopy(tu_mag)
            upr_node=deepcopy(nd)
        end
        tc_mag=lof_pnt2pnt_dist(nd.xy,cnt_xy)
        if (tc_mag<cnt_node_mag)
            cnt_node_mag=deepcopy(tc_mag)
            cnt_node=deepcopy(nd)
        end
    end


    #each owpp path
    lwr_paths=node[]
    upr_paths=node[]
    cnt_paths=node[]
    lwr_lths=0
    upr_lths=0
    cnt_lths=0
    power_sum=0
    wind_sum=wind()
    wind_sum.pu=zeros(Float64,length(buses[1].wnd.pu))
    wind_sum.ce=zeros(Float64,length(buses[1].wnd.ce))
    wind_sum.delta=0
    wind_sum.lf=0

    push!(upr_paths,deepcopy(as_Astar(pcc.node,upr_node,ocn.discretedom.nodes)))
    push!(lwr_paths,deepcopy(as_Astar(pcc.node,lwr_node,ocn.discretedom.nodes)))
    push!(cnt_paths,deepcopy(as_Astar(pcc.node,cnt_node,ocn.discretedom.nodes)))
    lwr_lths=lwr_lths+deepcopy(lwr_paths[length(lwr_paths)].G_cost)
    upr_lths=upr_lths+deepcopy(upr_paths[length(upr_paths)].G_cost)
    cnt_lths=cnt_lths+deepcopy(cnt_paths[length(cnt_paths)].G_cost)
    for opp in buses
        push!(upr_paths,deepcopy(as_Astar(opp.node,upr_node,ocn.discretedom.nodes)))
        push!(lwr_paths,deepcopy(as_Astar(opp.node,lwr_node,ocn.discretedom.nodes)))
        push!(cnt_paths,deepcopy(as_Astar(opp.node,cnt_node,ocn.discretedom.nodes)))
        lwr_lths=lwr_lths+deepcopy(lwr_paths[length(lwr_paths)].G_cost)
        upr_lths=upr_lths+deepcopy(upr_paths[length(upr_paths)].G_cost)
        cnt_lths=cnt_lths+deepcopy(cnt_paths[length(cnt_paths)].G_cost)
        wind_sum.pu=(wind_sum.pu.+deepcopy(opp.wnd.pu))
        wind_sum.ce=(wind_sum.ce.+deepcopy(opp.wnd.ce))
        wind_sum.delta=(wind_sum.delta+deepcopy(opp.wnd.delta))
        wind_sum.lf=(wind_sum.lf+deepcopy(opp.wnd.lf))
        power_sum=power_sum+deepcopy(opp.mva)
    end
    wind_sum.pu=(wind_sum.pu)./length(buses)
    wind_sum.ce=(wind_sum.ce)./length(buses)
    wind_sum.delta=(wind_sum.delta)/length(buses)
    wind_sum.lf=(wind_sum.lf)/length(buses)

    best_index=findmin([lwr_lths,cnt_lths,upr_lths])[2]
    paths=[lwr_paths,cnt_paths,upr_paths]
    nds=[lwr_node,cnt_node,upr_node]
    best_path=paths[best_index]
    best_node=paths[best_index][1]

    ojama_paths=node[]
    ojama_xys=xy[]
    for nd in ocn.owpps[owp.num+1:buses[length(buses)].num]
        dummy_path=as_Astar(nd,best_node,ocn.discretedom.nodes)
        push!(ojama_paths,deepcopy(dummy_path))
        goal=deepcopy(dummy_path.goal)
        while (dummy_path.parent.num != goal)
            dummy_path=dummy_path.parent
        end
        push!(ojama_xys,deepcopy(dummy_path.xy))
    end
    #adjust best connection point
    path2adjust=deepcopy(best_node)
    if (ymajorOWPP==true)
        room2adjust=deepcopy(best_node.xy.x)
        if (best_node.xy.x<=owp.node.xy.x)
            goal=deepcopy(path2adjust.goal)
            while path2adjust.num != goal
                path2adjust=path2adjust.parent
                if (path2adjust.xy.x<room2adjust)
                    room2adjust=deepcopy(path2adjust.xy.x)
                end
            end
            for xy in ojama_xys
                if (xy.x<best_node.xy.x && xy.x>=room2adjust)
                    best_node.xy.x=deepcopy(best_node.xy.x)
                end
            end
        else
            goal=deepcopy(path2adjust.goal)
            while path2adjust.num != goal
                path2adjust=path2adjust.parent
                if (path2adjust.xy.x>room2adjust)
                    room2adjust=deepcopy(path2adjust.xy.x)
                end
            end
            for xy in ojama_xys
                if (xy.x>best_node.xy.x && xy.x<=room2adjust)
                    best_node.xy.x=deepcopy(best_node.xy.x)
                end
            end
        end
    else
        room2adjust=deepcopy(best_node.xy.y)
        goal=deepcopy(path2adjust.goal)
        if (best_node.xy.y<=owp.node.xy.y)
            while path2adjust.num != goal
                path2adjust=path2adjust.parent
                if (path2adjust.xy.y<room2adjust)
                    room2adjust=deepcopy(path2adjust.xy.y)
                end
            end
            for xy in ojama_xys
                if (xy.y<best_node.xy.y && xy.y>=room2adjust)
                    best_node.xy.y=deepcopy(best_node.xy.y)
                end
            end
        else
            while path2adjust.num != goal
                path2adjust=path2adjust.parent
                if (path2adjust.xy.y>room2adjust)
                    room2adjust=deepcopy(path2adjust.xy.y)
                end
            end
            for xy in ojama_xys
                if (xy.y>best_node.xy.y && xy.y<=room2adjust)
                    best_node.xy.y=deepcopy(best_node.xy.y)
                end
            end
        end
    end
    domain_oss=deepcopy(ocn.discretedom.nodes)
    edges_oss=deepcopy(ocn.discretedom.edges)
    bs,best_node=opt_ifIn2Out(best_node,ocn)
    if bs==true
        best_node.num=deepcopy(length(domain_oss)+1)
        push!(domain_oss,best_node)
        lof_edgeifyOss(best_node,ocn,domain_oss,edges_oss)
    end
    return best_node,domain_oss,edges_oss,power_sum,wind_sum
end

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
        for xf in oss.xfmrs
            circ.cost=circ.cost+xf.costs.ttl
        end
    end
    for oss in circ.osss_mog
        for xf in oss.xfmrs
            circ.cost=circ.cost+xf.costs.ttl
        end
    end
    for xf in circ.pcc.xfmrs
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
