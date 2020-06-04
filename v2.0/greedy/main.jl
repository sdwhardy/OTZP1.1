include("functions.jl")
ocean2=deepcopy(ocn)
ocn=deepcopy(ocean2)
greedy_search(ocn)
i=7
beta=D[1]
function greedy_search(ocn)
    #make a copy of Tb_hv and Tb_mv into one set
    circuits_set=ocn.hv_circuits
    begin_at=0
    for (i,c) in enumerate(ocn.mv_circuits)
        if (begin_at==0 && c[1].decimal>6)
            begin_at=copy(i)
        end
        if (circuits_set[i][1].id != c[1].id)
            push!(circuits_set[i],ocn.mv_circuits[i][1])
        end
    end

    for i=begin_at:1:length(circuits_set)
        #decimal of index j
        j=circuits_set[i][1].decimal
        println("rank: "*string(j))
        #binary of index j
        j_bn=circuits_set[i][1].binary
        #Set D
        D=findall(x->x==1,j_bn)
        #proceed if enough OWPPs are included for breakdown
        if (length(D)>2)
            #copy original hv topology
            hv_orig=deepcopy(ocn.hv_circuits[i][1])
            #copy original mv topology
            mv_orig=deepcopy(ocn.mv_circuits[i][1])
            #index beta
            for beta in D
                beta_n_below,above_beta=split_j_at_beta(j_bn,beta[2])
                beta2gamma=findall(x->x==1,above_beta)#number of OWPPs above beta
                #proceed if more than 2
                if (length(beta2gamma)>1)
                    forwardHV,hv_remMV=copy_below_beta_equipment(deepcopy(hv_orig),beta_n_below,ocn,beta[2])
                    forwardMV,mv_remMV=copy_below_beta_equipment(deepcopy(mv_orig),beta_n_below,ocn,beta[2])
                    #Adjust MV OSS position if MV cable removed
                    if (mv_remMV==true)
                        forwardMV=add_temporary_hv_connections(forwardMV,ocn,beta[2])
                        forwardMV=finalize_circuit_layout(forwardMV,ocn,ks)
                        plot_circuit(forwardMV)
                        forwardMV,mv_remMV=copy_below_beta_equipment(forwardMV,beta_n_below,ocn,beta[2])
                        identical=false
                    else#if MV cable not removed check if HV and MV base are identical
                        #NOTE stopped here above works but needds to be tested with higher numbers of HV/MV/...cables
                        identical=check_pair(forwardHV,forwardMV)
                    end
                    #decompose circuits
                    if (identical==false)
                    #if (length(circuits_set[i])>1)
                        forwardHVTemp=[]
                        forwardMVTemp=[]
                        forwardHVTemp=opt_decomposeLO_solo(forwardHV, i_min_sum_bn, deepcopy(circuits_set),ocn,nsb_sum_bn)
                        forwardMVTemp=opt_decomposeLO_soloMV(forwardMV, i_min_sum_bn, deepcopy(circuits_set),ocn,nsb_sum_bn)
                        #if ((length(forwardHVTemp))>ceil(Int,opt_set()))
                        forwardHVTemp=opt_unicOnly(forwardHVTemp)
                            #forwardHVTemp=forwardHVTemp[1:ceil(Int,opt_set())]
                        #else
                        #end
                        #if ((length(forwardMVTemp))>ceil(Int,opt_set()))
                        forwardMVTemp=opt_unicOnly(forwardMVTemp)
                            #forwardMVTemp=forwardMVTemp[1:ceil(Int,opt_set())]
                        #else
                        #end
                        rnk=round(Int,forwardHV.decimal)
                        for ch in forwardHVTemp
                            insert_and_dedup!(circuits_set[i],deepcopy(ch))
                        end
                        for cm in forwardMVTemp
                            insert_and_dedup!(circuits_set[i],deepcopy(cm))
                        end
                        #circuits_set[round(Int,forwardHV.decimal)]=opt_unicOnly(circuits_set[round(Int,forwardHV.decimal)])
                        #circuits_set[round(Int,forwardHV.decimal)]=circuits_set[round(Int,forwardHV.decimal)][1:1]
                        #println("circuits_set: "*string(length(circuits_set[round(Int,forwardHV.decimal)])))

                        #circuits_set[i]=keep_unic(circuits_set[i])
                    else
                        forwardTemp=opt_decomposeLO_solo(forwardHV, i_min_sum_bn, deepcopy(circuits_set),ocn,nsb_sum_bn)
                        forwardTemp=opt_unicOnly(forwardTemp)
                        for cm in forwardTemp
                            insert_and_dedup!(circuits_set[i],deepcopy(cm))
                        end
                    end

                end
            end
        end
    end
    return begin_at, circuits_set
end
