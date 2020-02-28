
function opt_reAdjust_oss_1by1(systm,mv_rng,mv_cl)

    oss_x=Float64[]
    oss_y=Float64[]
    mog_x=Float64[]
    mog_y=Float64[]
    eps=10^-8
    for (ossORmog,junctions) in enumerate([systm.osss_owp,systm.osss_mog])
        for junction in junctions
            opt_bool=true
            for o2o_c in systm.oss2oss_cbls
                if (o2o_c.pth[length(o2o_c.pth)].num==junction.node.num)
                    opt_bool=false
                end
            end
            if (opt_bool=true)
                #create model and get OSS/MOG points
                m = Model(with_optimizer(Ipopt.Optimizer, print_level=1))
                #m = Model(solver=NLoptSolver(algorithm=:LD_MMA))
                #m = Model(with_optimizer(NLoptSolver,algorithm=LD_MMA))
                @variable(m,x,start=junction.node.xy.x)
                @variable(m,y,start=junction.node.xy.y)

                connections=Array{Tuple{Any,Any,Any,Any,Float64},1}()
                for mv_c in systm.owp_MVcbls
                    if (junction.node.num==mv_c.pth[length(mv_c.pth)].num)
                        push!(connections,(mv_c.pth[1].xy.x, mv_c.pth[1].xy.y, x, y, 2.52))
                        @constraint(m,((connections[length(connections)][1]-connections[length(connections)][3])^2+(connections[length(connections)][2]-connections[length(connections)][4])^2) <= mv_rng^2)
                    end
                end
                for hv_c in systm.owp_HVcbls
                    if (junction.node.num==hv_c.pth[1].num)
                        push!(connections,(x, y, hv_c.pth[length(hv_c.pth)].xy.x, hv_c.pth[length(hv_c.pth)].xy.y, hv_c.costs.ttl/hv_c.length))
                    elseif (junction.node.num==hv_c.pth[length(hv_c.pth)].num)
                        push!(connections,(hv_c.pth[1].xy.x, hv_c.pth[1].xy.y, x, y, hv_c.costs.ttl/hv_c.length))
                    end
                end
                for o2o_c in systm.oss2oss_cbls
                    if (junction.node.num==o2o_c.pth[1].num)
                        push!(connections,(x, y, o2o_c.pth[length(o2o_c.pth)].xy.x, o2o_c.pth[length(o2o_c.pth)].xy.y, o2o_c.costs.ttl/o2o_c.length))
                    elseif (junction.node.num==o2o_c.pth[length(o2o_c.pth)].num)
                        push!(connections,(o2o_c.pth[1].xy.x, o2o_c.pth[1].xy.y, x, y, o2o_c.costs.ttl/o2o_c.length))
                    end
                end
                for pcc_c in systm.pcc_cbls
                    if (junction.node.num==pcc_c.pth[1].num)
                        push!(connections,(x, y, pcc_c.pth[length(pcc_c.pth)].xy.x, pcc_c.pth[length(pcc_c.pth)].xy.y, pcc_c.costs.ttl/pcc_c.length))
                    end
                end
                #=for con in connections
                    println(con)
                end=#
                @NLobjective(m, Min,sum(((((connections[i][1]-connections[i][3])^2+(connections[i][2]-connections[i][4])^2)+eps)^(0.5))*connections[i][5] for i in 1:length(connections)))
                optimize!(m)
                junction.node.xy.x,junction.node.xy.y=(JuMP.value.(x),JuMP.value.(y))
                for mv_c in systm.owp_MVcbls
                    if (junction.node.num==mv_c.pth[length(mv_c.pth)].num)
                        mv_c.pth[length(mv_c.pth)].xy.x=deepcopy(junction.node.xy.x)
                        mv_c.pth[length(mv_c.pth)].xy.y=deepcopy(junction.node.xy.y)
                    end
                end
                for hv_c in systm.owp_HVcbls
                    if (junction.node.num==hv_c.pth[1].num)
                        hv_c.pth[1].xy.x=deepcopy(junction.node.xy.x)
                        hv_c.pth[1].xy.y=deepcopy(junction.node.xy.y)
                    elseif (junction.node.num==hv_c.pth[length(hv_c.pth)].num)
                        hv_c.pth[length(hv_c.pth)].xy.x=deepcopy(junction.node.xy.x)
                        hv_c.pth[length(hv_c.pth)].xy.y=deepcopy(junction.node.xy.y)
                    end
                end
                for o2o_c in systm.oss2oss_cbls
                    if (junction.node.num==o2o_c.pth[1].num)
                        o2o_c.pth[1].xy.x=deepcopy(junction.node.xy.x)
                        o2o_c.pth[1].xy.y=deepcopy(junction.node.xy.y)
                    elseif (junction.node.num==o2o_c.pth[length(o2o_c.pth)].num)
                        o2o_c.pth[length(o2o_c.pth)].xy.x=deepcopy(junction.node.xy.x)
                        o2o_c.pth[length(o2o_c.pth)].xy.y=deepcopy(junction.node.xy.y)
                    end
                end
                for pcc_c in systm.pcc_cbls
                    if (junction.node.num==pcc_c.pth[1].num)
                        pcc_c.pth[1].xy.x=deepcopy(junction.node.xy.x)
                        pcc_c.pth[1].xy.y=deepcopy(junction.node.xy.y)
                    end
                end
                if (ossORmog==1)
                    push!(oss_x,deepcopy(junction.node.xy.x))
                    push!(oss_y,deepcopy(junction.node.xy.y))
                else
                    push!(mog_x,deepcopy(junction.node.xy.x))
                    push!(mog_y,deepcopy(junction.node.xy.y))
                end
            else
            end
        end
    end
    return [mog_x,mog_y,oss_x,oss_y]
end



@time bsf, ocn=testing_loop()

ppf_equipment_OSS_MOG(ocn,bsf[1])
ppf_equipment_OSS_MOG(ocean_bsf,savC)


######################################## depricated
#=
S=500
kv=220
l=10
sze=500
nm=2.
cst=100
wp=ocean.owpps[1].wnd
#ks=ocean.finance
cst=50.30172
cstF_hVrng_o2o(l,S,wp,ocean.finance,ocean.sys,sze,kv,nm,cst)#NOT WORKING AT ALLLLLL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function cstF_hVrng_o2o(l,S,wp,ks,sys,sze,kv,nm,cst)
    hv0,hv1=cst_nextSizeUp(l,S,sze,kv,nm)
    cst_hv0=1
    cst_hv1=Inf
    cst_km=cst/l+(cst/l)*0.005
    cst_hv0_km=cst_hv0/l
    l_lim=deepcopy(2*l)
    while ((cst_hv0 < cst_hv1) && (cst_km>cst_hv0_km) && (l<l_lim))
        l=l+1
        cst_hv0=cst_setCable(l,S,wp,ks,hv0)
        cst_hv1=cst_setCable(l,S,wp,ks,hv1)
        cst_hv0_km=cst_hv0/l
        println(string(cst_hv0)*" - "*string(cst_hv1))
    end
    return l-1
end

function cst_setCable(l,S,wp,ks,value)
    value.costs.cbc=cstF_cbl_cpx(value)#capex
    value.costs.qc=cstF_cbl_qo2o(value,ks)#cost of compensastion
    value.costs.rlc=cstF_acCbl_rlc(value,S,ks,wp)#cost of losses
    value.costs.cm=cstF_eqp_cm(value,ks)#corrective maintenance
    value.costs.eens=eensF_eqp_eens(value,S,ks,wp)#eens calculation
    value.costs.ttl=cstF_cbl_sum(value.costs)#totals the cable cost
end

function wnd_setPandWonOss(ocn,circs)
    for (j,circ) in enumerate(circs)
        for mv_cbl in circ.owp_MVcbls
            keta_os=mv_cbl.pth[length(mv_cbl.pth)].num
            keta_owp=mv_cbl.pth[1].num
            for (i,os) in enumerate(circ.osss_owp)
                if (os.node.num==keta_os)
                    println(string(j)*" changed")
                    circ.osss_owp[i].mva=ocn.owpps[keta_owp].mva
                    circ.osss_owp[i].wnd=ocn.owpps[keta_owp].wnd
                end
            end
        end
    end
    return circs
end

=#


#=circs=bsf_mv
circ=circs[7]
hv_cbl=circ.owp_HVcbls[1]
mv_cbl=circ.owp_MVcbls[1]
os=circ.osss_owp[1]
function opt_cableLimits(circs,ocn)
    cableMX_dict = Dict()
    for circ in circs
        for hv_cbl in circ.owp_HVcbls
            S=0
            wp=wind()
            for os in circ.osss_owp
                if (os.node.num==hv_cbl.pth[1].num)
                    S=os.mva
                    wp=os.wnd
                end
            end
            if (S==0)
                S=circ.owpps[1].mva
                wp=circ.owpps[1].wnd
            end
            hv_name=string(S)*"_"*string(round(Int32,hv_cbl.length))*"_"*string(hv_cbl.elec.volt)*"_"*string(hv_cbl.size)*"_"*string(hv_cbl.num)
            if (haskey(dict, hv_name))
                l_mx=cstF_hVrng_o2o(hv_cbl.length,S,wp,ocn.finance,ocn.sys,hv_cbl.size,hv_cbl.elec.volt,hv_cbl.num,hv_cbl.costs.ttl)
                dict[hv_name]=l_mx
            end
        end
        for o2o_cbl in circ.oss2oss_cbls
            S=0
            wp=wind()
            for mog in circ.osss_mog
                if (mog.node.num==o2o_cbl.pth[1].num)
                    S=mog.mva
                    wp=mog.wnd
                end
            end
            if (S==0)
                S=circ.owpps[1].mva
                wp=circ.owpps[1].wnd
            end
            o2o_name=string(S)*"_"*string(round(Int32,o2o_cbl.length))*"_"*string(o2o_cbl.elec.volt)*"_"*string(o2o_cbl.size)*"_"*string(o2o_cbl.num)
            if (haskey(dict, o2o_name))
                l_mx=cstF_hVrng_o2o(o2o_cbl.length,S,wp,ocn.finance,ocn.sys,o2o_cbl.size,o2o_cbl.elec.volt,o2o_cbl.num,o2o_cbl.costs.ttl)
                dict[hv_name]=l_mx
            end
        end
        for pcc_cbl in circ.pcc_cbls
            S=0
            wp=wind()
            for mog in circ.osss_mog
                if (mog.node.num==pcc_cbl.pth[1].num)
                    S=mog.mva
                    wp=mog.wnd
                end
            end
            if (S==0)
                S=circ.owpps[1].mva
                wp=circ.owpps[1].wnd
            end
            pcc_name=string(S)*"_"*string(round(Int32,pcc_cbl.length))*"_"*string(pcc_cbl.elec.volt)*"_"*string(pcc_cbl.size)*"_"*string(pcc_cbl.num)
            if (haskey(dict, pcc_name))
                l_mx=cstF_hVrng_o2p(pcc_cbl.length,S,wp,ocn.finance,ocn.sys,pcc_cbl.size,pcc_cbl.elec.volt,pcc_cbl.num)
                dict[pcc_name]=l_mx
            end
        end
    end
    return dict
end=#
