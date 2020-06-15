include("functions.jl")
#=ocn=deepcopy(ocean)
i=7
beta=D[1]
r=greedy_search(ocn)
original=originals[1]
circuits_set=deepcopy(TH)
original=originals[1]
tbeta=tbetas[1]=#
function greedy_search(ocn)
    circuits_set, begin_at = pre_filter_circuits(ocn.hv_circuits,ocn.mv_circuits)
    ks=get_Cost_Data()

    for i=begin_at:1:length(circuits_set)
        #decimal of index j
        j=circuits_set[i][1].decimal
        println(" owpps: "*string(circuits_set[i][1].binary)*" ("*string(j)*")")
        #binary of index j
        j_bn=circuits_set[i][1].binary
        #Set D
        D=findall(x->x==1,j_bn)
        #proceed if enough OWPPs are included for breakdown
        if (length(D)>2)
            #copy originals hv topology
            originals=circuit[]
            for circ in circuits_set[i]
                push!(originals,deepcopy(circ))
            end
            #index beta
            TH_jbeta=circuits_set[i]
            println("old circ "*string(circuits_set[i][1].decimal)*" costs: "*string(circuits_set[i][1].cost))
            for beta in D
                beta_n_below,above_beta=split_j_at_beta(j_bn,beta[2])
                beta2gamma=findall(x->x==1,above_beta)#number of OWPPs above beta
                #proceed if more than 2
                tbetas=circuit[]
                if (length(beta2gamma)>1)
                    for original in originals
                        tbeta_temp,remMV=copy_below_beta_equipment(deepcopy(original),beta_n_below,ocn,beta[2])
                    #Adjust MV OSS position if MV cable removed
                        if (remMV==true)
                            tbeta_temp=add_temporary_hv_connections(tbeta_temp,ocn,beta[2],ks)
                            tbeta_temp=finalize_circuit_layout(tbeta_temp,ocn,ks)
                            tbeta_temp,remMV=copy_below_beta_equipment(tbeta_temp,beta_n_below,ocn,beta[2])
                        end
                        if (length(tbetas)>1)
                            if (check_if_identical(tbetas[1],tbeta_temp)==false)
                                push!(tbetas,tbeta_temp)
                            end
                        else
                            push!(tbetas,tbeta_temp)
                        end
                    end

                    #decompose circuits
                    tH_mvhvs=Vector{Array{circuit,1}}()
                    for (tbeta_1,tbeta) in enumerate(tbetas)
                        tH_mvhv=find_tH(tbeta, beta_n_below, above_beta, deepcopy(circuits_set[1:i-1]),ocn,ks)
                        push!(tH_mvhvs,tH_mvhv)
                    end

                    for tH_mvhv in tH_mvhvs
                        for tH in tH_mvhv
                            if (is_a_unique_entry(tH, TH_jbeta))
                                insert_and_dedup!(TH_jbeta, tH)
                                #TH_jbeta=adjust_length_of_q(TH_jbeta,length_of_THs_i())
                            end
                        end
                    end
                end
            end
            if (i<length(circuits_set))
                println("length_b4: "*string(length(TH_jbeta)))
                TH_jbeta=check_circuits_2_keep(TH_jbeta,circuits_set[i+1:length(circuits_set)],ocn)
                println("length_afta: "*string(length(TH_jbeta)))
            end
            circuits_set[i]=TH_jbeta
            println("new circ "*string(circuits_set[i][1].decimal)*" costs: "*string(circuits_set[i][1].cost))
        end
    end
    return circuits_set
end

function check_circuits_2_keep(TH_jbeta,circuits_set,ocn)
    TH_jbeta2_keep=TH_jbeta[1:1]
    base_aft_mog=TH_jbeta[1].mog[1].node
    base_cost=TH_jbeta[1].cost
    base_1s=findall(x->x==1,TH_jbeta[1].binary)
    forward_circs=Vector{Array{circuit,1}}()
    for circuit in circuits_set
        circ_1s=findall(x->x==1,circuit[1].binary)
        if issubset(base_1s,circ_1s)
            push!(forward_circs,circuit)
        end
    end

    for th_jb in TH_jbeta[2:length(TH_jbeta)]
        keep=false
        cost_dif=th_jb.cost-base_cost
        for mog_circs in forward_circs
            for mog_circ in mog_circs
                km_dif=euclidian_distance(mog_circ.mog[length(mog_circ.mog)].node.xy,base_aft_mog.xy)-euclidian_distance(mog_circ.mog[length(mog_circ.mog)].node.xy,th_jb.mog[length(th_jb.mog)].node.xy)
                if ((cost_dif-km_dif*th_jb.PCCcbls[length(th_jb.PCCcbls)].costs.perkm_ttl)<=0)
                    keep=true
                    @goto shes_a_keeper_son
                end
            end
        end
        if (keep==true)
            @label shes_a_keeper_son
            push!(TH_jbeta2_keep,th_jb)
        end
    end

    return TH_jbeta2_keep
end


#=save("v2.0/tempfiles/ocean/tH_hv.jld2", "tH_hv",tH_hv)
save("v2.0/tempfiles/ocean/tH_mv.jld2", "tH_mv",tH_mv)

for (i,t) in enumerate(tH_mv)
    println(string(i)*" - "*string(t.cost))
end
plotly()
p=plot()
plot_circuit(p,tH_mv[5])
gui()=#
