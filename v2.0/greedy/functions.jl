#sets size of queues
function length_of_TEs_i()
    return 1
end

function length_of_TEs_fr()
    return 1
end

function length_of_THs_i()
    return 1
end

function length_of_THs_fr()
    return 1000
end

function adjust_length_of_q(q,l)
    if (length(q)>l)
        q=q[1:l]
    else
    end
    return q
end

#keeps priority queues properly ordered
insert_and_dedup!(v::Vector, x) = (splice!(v, searchsorted(v,x,by = v -> v.cost), [x]); v)

insert_and_dedup_dict!(v::Vector, x) = (splice!(v, searchsorted(v,x,by = v -> v["cost"]), [x]); v)

function pre_filter_circuits(hv_circuits,mv_circuits)
    circuit_set=hv_circuits
    begin_at=0
    #keep all HV circuits that have HV cables longer than the MV connections
    for (hv_i,hv_circ) in enumerate(hv_circuits)
        if (begin_at==0 && hv_circ[1].decimal>6)
            begin_at=copy(hv_i)
        end
        owpp_oss_pairs=Array{Tuple,1}()
        for mv_c in hv_circ[1].MVcbls
            for owp in hv_circ[1].owpps
                if (owp.node.num==mv_c.path[1].num)
                    push!(owpp_oss_pairs,(owp,mv_c.path[2].num))
                end
            end
        end

        for hv_c in hv_circ[1].HVcbls
            for owpp_oss_pair in owpp_oss_pairs
                if (owpp_oss_pair[2]==hv_c.path[1].num)
                    if (owpp_oss_pair[1].mv_zone>hv_c.length)
                        @goto remove_the_hv_circuit
                    end
                end
            end
        end
        if (false)
            @label remove_the_hv_circuit
            circuit_set[hv_i]=circuit[]
        end
    end
    #keep MV circuits if a direct MV connection exists
    for (mv_i,mv_circ) in enumerate(mv_circuits)
        if ((length(mv_circ[1].MVcbls)>length(mv_circ[1].HVcbls)) && (length(mv_circ[1].MVcbls)>1))
            push!(circuit_set[mv_i],mv_circ[1])
        end
    end
    return circuit_set, begin_at
end

function is_a_valid_combination(bn0,bn1,l)
    cmb=true
    bn01=zeros(Int8,l)
    for one0 in bn0
        for one1 in bn1
            if (one0==one1)
                cmb=false
            end
        end
    end
    if (cmb==true)
        for one0 in bn0
            bn01[one0]=1
        end
        for one1 in bn1
            bn01[one1]=1
        end
    end
    return cmb, bn01
end

function is_a_unique_entry(new_circ, circs)
    uniq=true
    if (length(circs)>0)
        nearest_circs=circuit[]
        for circ in circs
            if (circ.cost<=new_circ.cost)
                nearest_circs=[circ]
            end
            if (circ.cost>new_circ.cost)
                push!(nearest_circs,circ)
                break
            end
        end
        for circ in nearest_circs
            uniq=same_circuit_cost_and_id_test(new_circ,circ)
        end
    end
    return uniq
end

function is_a_unique_entry_dict(new_circ, circs)
    uniq=true
    if (length(circs)>0)
        nearest_circs=Vector{Dict{String,Any}}()
        for circ in circs
            if (circ["cost"]<=new_circ["cost"])
                nearest_circs=[new_circ]
            end
            if (circ["cost"]>new_circ["cost"])
                push!(nearest_circs,circ)
                break
            end
        end
        for circ in nearest_circs
            uniq=same_circuit_cost_and_id_test_dict(new_circ,circ)
        end
    end
    return uniq
end

function same_circuit_cost_and_id_test(new_circ,circ)
    uniq=true
    err=10e-2
    if (new_circ.cost-err<circ.cost && new_circ.cost+err>circ.cost)
        uniq=same_circuit_id_test(new_circ,circ)
    end
    return uniq
end

function same_circuit_cost_and_id_test_dict(new_circ,circ)
    uniq=true
    err=10e-2
    if (new_circ["cost"]-err<circ["cost"] && new_circ["cost"]+err>circ["cost"])
        uniq=same_circuit_id_test_dict(new_circ,circ)
    end
    return uniq
end

function same_circuit_id_test_dict(new_circ,circ)
    uniq=true
    str0=deepcopy(new_circ["id"])
    str1=deepcopy(circ["id"])
    tst0=split(str0, "_")
    tst1=split(str1, "_")
    a=sort!(tst0)
    b=sort!(tst1)
    if (a==b)
        uniq=false
    end
    return uniq
end

function same_circuit_id_test(new_circ,circ)
    uniq=true
    str0=deepcopy(new_circ.id)
    str1=deepcopy(circ.id)
    tst0=split(str0, "_")
    tst1=split(str1, "_")
    a=sort!(tst0)
    b=sort!(tst1)
    if (a==b)
        uniq=false
    end
    return uniq
end

function split_j_at_beta(j,beta)
    beta_n_below=Int8[]
    above_beta=Int8[]
    for (position,bit) in enumerate(j)
        if (position <= beta)
            push!(beta_n_below,bit)
            push!(above_beta,0)
        elseif (position > beta)
            push!(beta_n_below,0)
            push!(above_beta,bit)
        end
    end
    return beta_n_below,above_beta
end

function convert_bin2dec(bn)
    dec=0.0
    for (i,bt) in enumerate(bn)
        dec=dec+bt*(2^(i-1))
    end
    return dec
end

function convert_dec2bin(dec)
    bn=Int8[]
    while (dec != 0)
        push!(bn,mod(dec,2))
        dec=floor(Int32, dec/2)
    end
    return bn
end

function copy_below_beta_equipment(crc,beta_n_below,ocn,beta)
    removedMVmain=false
    oss_system=circuit()
    oss_system.binary=crc.binary
    oss_system.pcc=crc.pcc
    oss_system.decimal=crc.decimal
    oss_system.owpps=crc.owpps
    oss_system.base=crc.base
    oss_system.wnd=crc.wnd
    oss_system.mva=crc.mva
    oss_system.id=crc.id

    #Take base OWP,PCC connections and MOG without transformers
    push!(oss_system.mog,crc.mog[1])
    oss_system.mog[1].xfmrs=transformer[]
    push!(oss_system.PCCcbls,crc.PCCcbls[1])
    mog=oss_system.mog[1].node.num
    oss_mog=Int32[]
    for (owp_i,bn) in enumerate(beta_n_below)
        if (bn==1)
            owp=ocn.owpps[owp_i].node.num
            for mv_c in crc.MVcbls
                if (mv_c.path[1].num==owp)
                    push!(oss_system.MVcbls,mv_c)
                    if (mv_c.path[length(mv_c.path)].num!=mog)
                        #copy numbers of nodes at end of MV that aren't the terminous MOG
                        push!(oss_mog,deepcopy(mv_c.path[length(mv_c.path)].num))
                    end
                end
                if (mv_c.path[1].num>beta && mv_c.path[length(mv_c.path)].num==mog)
                    removedMVmain=true
                end
            end
        end
    end

    for hv_c in crc.HVcbls
        for oss_num in oss_mog
            if (hv_c.path[1].num==oss_num)
                push!(oss_system.HVcbls,hv_c)
                if (hv_c.path[length(hv_c.path)].num != mog)
                    push!(oss_mog,deepcopy(hv_c.path[length(hv_c.path)].num))
                end
            end
        end
    end

    for mvhv_oss in crc.oss
        for oss_num in oss_mog
            if (mvhv_oss.node.num==oss_num)
                push!(oss_system.oss,mvhv_oss)
            end
        end
    end

    for mvhv_mog in crc.mog[2:length(crc.mog)]
        for oss_num in oss_mog
            if (mvhv_mog.node.num==oss_num)
                push!(oss_system.mog,mvhv_mog)
            end
        end
    end
    for o2o_c in crc.O2Ocbls
        for oss_num in oss_mog
            if (o2o_c.path[1].num==oss_num)
                push!(oss_system.O2Ocbls,o2o_c)
            end
        end
    end
    total_circuit_cost(oss_system)
    return oss_system,removedMVmain
end

function add_temporary_hv_connections(forwardMV,ocn,beta,ks)
    destination=forwardMV.mog[1].node
    fake_oss=0
    for owp in forwardMV.owpps
        if (owp.node.num>beta)
            wnd_power=create_wind_and_power_dict()
            km=euclidian_distance(owp.node.xy,destination.xy)
            owpp2oss_connections=owpp_acdc_to_oss(km,owp,destination,ocn.database,ks,true)
            cheapest=argmin([owpp2oss_connections["300.0"]["cost"],owpp2oss_connections["220.0"]["cost"],owpp2oss_connections["400.0"]["cost"]])
            hv_connection=[owpp2oss_connections["300.0"],owpp2oss_connections["220.0"],owpp2oss_connections["400.0"]][cheapest]
            forwardMV,wnd_power=best_owpp2oss(hv_connection,forwardMV,owp,ocn.database,ks,wnd_power,ocn)

            fake_oss=fake_oss+1
        end
    end
    forwardMV=total_circuit_cost(forwardMV)
    ocn.num=ocn.num-fake_oss
    return forwardMV
end


function check_if_identical(c0,c1)
    crcCopy=false
    if (length(c0.mog)==length(c1.mog) && length(c0.oss)==length(c1.oss) && length(c0.MVcbls)==length(c1.MVcbls) && length(c0.HVcbls)==length(c1.HVcbls) && length(c0.O2Ocbls)==length(c1.O2Ocbls) && length(c0.PCCcbls)==length(c1.PCCcbls))
        crcCopy=true
    end
    return crcCopy
end




function find_tH(tc, beta_n_below, above_beta, circuits_set, ocn, ks)
    tEs=find_best_tEs(above_beta, circuits_set, tc, ocn, ks)
    #println("tEs: "*string(length(tEs)))
    cs=circuit[]
    if (length(tEs)>0)
        for tE in tEs
            wnd_power=create_wind_and_power_dict()
            single_parent=deepcopy(tc)
            single_parent.id=single_parent.id*"_"*tE.id
            #copy the base wind and power values

            for mv_c in single_parent.MVcbls
                if (mv_c.path[2].num==single_parent.mog[1].node.num)
                    wnd_power=track_wind_and_power(wnd_power,mv_c)
                end
            end
            for hv_c in single_parent.HVcbls
                if (hv_c.path[2].num==single_parent.mog[1].node.num)
                    wnd_power=track_wind_and_power(wnd_power,hv_c)
                end
            end
            for hv_c in single_parent.O2Ocbls
                if (hv_c.path[2].num==single_parent.mog[1].node.num)
                    wnd_power=track_wind_and_power(wnd_power,hv_c)
                end
            end
            #copy the upstream wind and power values
            single_parent, wnd_power = transfer_2_circuits_equipment(single_parent,tE,wnd_power)
            if (single_parent.mog[1].kV==300.0)
                single_parent.mog[1].conv=set_oss_converters(single_parent.mog[1].kV,wnd_power,ocn.database,ks)
            else
            end
            #find transformers
            single_parent.mog[1].xfmrs=set_oss_transformers(single_parent.mog[1].kV,wnd_power,ocn.database,ks)
            #set platform structure
            single_parent=set_mog_platform(single_parent,ks)
            single_parent=total_circuit_cost(single_parent)
            if (is_a_unique_entry(single_parent, cs))
                insert_and_dedup!(cs,single_parent)
            end
        end
    end
    return cs
end


function find_best_tEs(above_beta, circuits_set, tc, ocn, ks)
    #println(string(top_bin2dec(aft_bn))*" :aft_bn: "*string(aft_bn))
    above_beta_1s=findall(x->x==1,above_beta)
    if (length(above_beta_1s)>1)
        #println("find_set_TE_parts: ")
        TE_parts, mid_point_index=find_set_TE_parts(above_beta, circuits_set, above_beta_1s, tc, ocn.owpps)
        #println("change_tEs_bases_pcc2oss: ")
        TE_parts=change_tEs_bases_pcc2oss(TE_parts, tc, ocn, ks)
        #println("best_tEs_from_TE_parts: ")
        best_tEs=best_tEs_from_TE_parts(TE_parts, above_beta, mid_point_index, ocn,ks)
    end
    return best_tEs
end

#NOTE just before going to exam made it to this point
function best_tEs_from_TE_parts(TE_parts, above_beta, mid_point_index, ocn,ks)
    #find all the complete systems
    complete_systems=circuit[]
    above_beta_dec=round(Int32,convert_bin2dec(above_beta))
    for TE_part in TE_parts[length(TE_parts)]
        if (TE_part.decimal==above_beta_dec)
            if (is_a_unique_entry(TE_part, complete_systems))
                insert_and_dedup!(complete_systems, TE_part)
                complete_systems=adjust_length_of_q(complete_systems,length_of_TEs_fr())
            end
        end
    end
    #roll up the components keeping the best and storing the complete systems
    for (TE_part_set_i,TE_part_set) in enumerate(TE_parts[1:mid_point_index])
        #find OWPPs included in first circuit
        TE_part_set_bin0,TE_part_set_bin1=split_j_at_beta(TE_part_set[length(TE_part_set)].binary,0)
        TE_part_1s=findall(x->x==1,TE_part_set_bin1)
        for (TE_part_i,TE_part) in enumerate(TE_part_set)
            for TE_part_next_set in TE_parts[TE_part_set_i+1:length(TE_parts)]
                #find OWPPs included in next circuit
                TE_part_next_set_bin0,TE_part_next_set_bin1=split_j_at_beta(TE_part_next_set[length(TE_part_next_set)].binary,0)
                TE_part_next_1s=findall(x->x==1,TE_part_next_set_bin1)
                #check if they are valid to combine
                valid_combo, new_binary=is_a_valid_combination(TE_part_1s,TE_part_next_1s,length(above_beta))
                if (valid_combo)
                    for (TE_part_next_i,TE_part_next) in enumerate(TE_part_next_set)
                        #find the index of the combined circuit
                        TE_parts_i=0
                        rnk=floor(Int64,convert_bin2dec(new_binary))
                        for (rnk_i,tep_set) in enumerate(TE_parts)
                            if (tep_set[1].decimal==rnk)
                                TE_parts_i=deepcopy(rnk_i)
                                break
                            end
                        end

                        if (TE_parts_i != 0)
                            #combine the 2 circuits
                            combined_system=combine_TE_parts(deepcopy(TE_parts[TE_parts_i][1]),deepcopy(TE_part),deepcopy(TE_part_next),ocn,ks)
                            if (above_beta_dec==rnk)#store the new full rank system
                                if (is_a_unique_entry(combined_system, complete_systems))
                                    insert_and_dedup!(complete_systems, combined_system)
                                    complete_systems=adjust_length_of_q(complete_systems,length_of_TEs_fr())
                                end
                            else#store the new partial rank system
                                if (is_a_unique_entry(combined_system, TE_parts[TE_parts_i]))
                                    insert_and_dedup!(TE_parts[TE_parts_i], combined_system)
                                    TE_parts[TE_parts_i]=adjust_length_of_q(TE_parts[TE_parts_i],length_of_TEs_i())
                                end
                            end
                        end
                    end
                end

            end
        end
    end
    return complete_systems
end


function FR_from_TH(TH_parts,FR_bin,ocn)
    ks=get_Cost_Data()
    gomi,above_beta=split_j_at_beta(FR_bin,0)
    above_beta_dec=round(Int32,convert_bin2dec(above_beta))
    mid_point_index=0
    for (TH_part_i,TH_part) in enumerate(TH_parts)
        if (TH_part[1].decimal<above_beta_dec/2+1)
            mid_point_index=deepcopy(TH_part_i)
        end
    end
    println(mid_point_index)
    #find all the complete systems
    complete_systems=circuit[]

    for TH_part in TH_parts[length(TH_parts)]
        if (TH_part.decimal==above_beta_dec)
            if (is_a_unique_entry(TH_part, complete_systems))
                insert_and_dedup!(complete_systems, TH_part)
                complete_systems=adjust_length_of_q(complete_systems,length_of_THs_fr())
            end
        end
    end

    #roll up the components keeping the best and storing the complete systems
    for (TH_part_set_i,TH_part_set) in enumerate(TH_parts[1:mid_point_index])
        #find OWPPs included in first circuit
        TH_part_set_bin0,TH_part_set_bin1=split_j_at_beta(TH_part_set[length(TH_part_set)].binary,0)
        TH_part_1s=findall(x->x==1,TH_part_set_bin1)
        for (TH_part_i,TH_part) in enumerate(TH_part_set)
            for TH_part_next_set in TH_parts[TH_part_set_i+1:length(TH_parts)]
                #find OWPPs included in next circuit
                TH_part_next_set_bin0,TH_part_next_set_bin1=split_j_at_beta(TH_part_next_set[length(TH_part_next_set)].binary,0)
                TH_part_next_1s=findall(x->x==1,TH_part_next_set_bin1)
                #check if they are valid to combine
                valid_combo, new_binary=is_a_valid_combination(TH_part_1s,TH_part_next_1s,length(above_beta))
                if (valid_combo)
                    for (TH_part_next_i,TH_part_next) in enumerate(TH_part_next_set)
                        #find the index of the combined circuit
                        rnk=floor(Int64,convert_bin2dec(new_binary))
                        combined_system=combine_TH_parts(deepcopy(TH_part),deepcopy(TH_part_next),reshape(new_binary, (1,length(new_binary))),rnk,ocn,ks)
                        if (above_beta_dec==rnk)#store the new full rank system
                            if (is_a_unique_entry(combined_system, complete_systems))
                                insert_and_dedup!(complete_systems, combined_system)
                                complete_systems=adjust_length_of_q(complete_systems,length_of_THs_fr())
                            end
                        else#store the new partial rank system
                            TH_parts_i=0#find the index of the combined circuit
                            TH_parts_spot=0
                            for (rnk_i,THp_set) in enumerate(TH_parts)
                                if (THp_set[1].decimal==rnk)
                                    TH_parts_i=deepcopy(rnk_i)
                                    break
                                elseif (THp_set[1].decimal>rnk)
                                    TH_parts_spot=deepcopy(rnk_i)
                                    break
                                end
                            end
                            if (TH_parts_i != 0)
                                if (is_a_unique_entry(combined_system, TH_parts[TH_parts_i]))
                                    insert_and_dedup!(TH_parts[TH_parts_i], combined_system)
                                    TH_parts[TH_parts_i]=adjust_length_of_q(TH_parts[TH_parts_i],length_of_THs_i())
                                end
                            elseif (TH_parts_spot != 0)
                                insert!(TH_parts,TH_parts_spot,[combined_system])
                            else
                                push!(TH_parts,[combined_system])
                            end
                        end
                    end
                end

            end
        end
    end
    return complete_systems
end

function fr_from_th(ths,ocn)
    ks=get_Cost_Data()
    fr_bin=ths[length(ths)]["bin"]
    fr_dec=ths[length(ths)]["dec"]
    mid_point_index=2^(length(fr_bin)-1)
    #find all the complete systems
    complete_systems=Vector{Dict{String,Any}}()
    insert_and_dedup!(complete_systems, ths[length(ths)])
    for (th0_i,th0) in enumerate(ths[1:mid_point_index])
        #find OWPPs included in first circuit
        th0_1s=findall(x->x==1,th0["bin"])
        for (th1_i,th1) in enumerate(ths[th0_i+1:length(ths)])
            #find OWPPs included in second circuit
            th1_1s=findall(x->x==1,th1["bin"])
            if (can_combine(th0_1s,th1_1s))
                th01=combine_ths(deepcopy(th0),deepcopy(th1),ocn,ks)
                println(fr_dec)
                println(th01["dec"])
                if (th01["dec"]==fr_dec)
                    if (is_a_unique_entry_dict(th01, complete_systems))
                        insert_and_dedup_dict!(complete_systems,deepcopy(th01))
                        #complete_systems=adjust_length_of_q(complete_systems,length_of_THs_fr())
                    end
                else
                    if (th01["cost"]<ths[round(Int64,th01["dec"])]["cost"])
                        ths[round(Int64,th01["dec"])]=deepcopy(th01)
                    end
                end
            end
        end
    end
    return complete_systems
end

function combine_ths(th0,th1,ocn,ks)
    th01=Dict{String,Any}()
    push!(th01,("mva"=>Float64[]))
    push!(th01,("dec"=>0.0))
    th01["bin"]=convert_dec2bin(2^length(ocn.owpps)-1)
    push!(th01,("cost"=>Inf))
    push!(th01,("xfo_cost"=>0.0))
    push!(th01,("conv_cost"=>0.0))
    push!(th01,("kV"=>Float64[]))
    push!(th01,("dec_set"=>Float64[]))
    push!(th01,("id"=>""))

    th01["dec"]=th0["dec"]+th1["dec"]
    th01["bin"]=convert_dec2bin(th01["dec"])
    #powers
    for mva in th0["mva"]
        push!(th01["mva"],mva)
    end
    for mva in th1["mva"]
        push!(th01["mva"],mva)
    end
    #kvs
    for kV in th0["kV"]
        push!(th01["kV"],kV)
    end
    for kV in th1["kV"]
        push!(th01["kV"],kV)
    end
    #dec_set
    for dec_set in th0["dec_set"]
        push!(th01["dec_set"],dec_set)
    end
    for dec_set in th1["dec_set"]
        push!(th01["dec_set"],dec_set)
    end

    #set new transformer and converter levels
    s220=0
    s400=0
    s300=0
    s220_cost=0
    s400_cost=0
    s300_costx=0
    s300_costc=0
    for (i,kV) in enumerate(th01["kV"])
        if (kV==220.0)
            s220=s220+th01["mva"][i]
        elseif (kV==400.0)
            s400=s400+th01["mva"][i]
        else
            s300=s300+th01["mva"][i]
        end
    end

    if (s220!=0.0 && ocn.pcc.kV!=220.0)
        xfo220=deepcopy(ocn.database["transformers"][string(s220)]["onshore"])
        s220_cost=xfo220.costs.ttl
    end
    if (s400!=0.0 && ocn.pcc.kV!=400.0)
        xfo400=deepcopy(ocn.database["transformers"][string(s400)]["onshore"])
        s400_cost=xfo400.costs.ttl
    end
    if (s300!=0.0)
        #converters
        conv300=deepcopy(ocean.database["converters"][string(s300)]["onshore"])
        s300_costc=conv300.costs.ttl
        #transformers
        xfo300=deepcopy(ocn.database["transformers"][string(s300)]["onshore"])
        s300_costx=xfo300.costs.ttl
    end


    th01["xfo_cost"]=s220_cost+s400_cost+s300_costx
    th01["conv_cost"]=s300_costc
    th01["cost"]=th0["cost"]+th1["cost"]-th0["xfo_cost"]-th1["xfo_cost"]-th0["conv_cost"]-th1["conv_cost"]+th01["xfo_cost"]+th01["conv_cost"]
    th01["id"]=th0["id"]*"_"*th1["id"]
    return th01
end

function can_combine(th0_1s,th1_1s)
    ok=false
    for th1_1 in th0_1s
        ok=issubset([th1_1], th1_1s)
        if (ok)
            break
        end
    end
    return !ok
end

function TH_reduction(TH,owpps)
    max_dec=2^owpps-1
    dict_of_connections=Array{Dict{String,Any},1}()
    counting_number=1
    for th in TH
        dict_entry=Dict{String,Any}()
        if (th[1].decimal==counting_number)
            push!(dict_entry,("mva"=>[th[1].PCCcbls[1].mva]))
            push!(dict_entry,("dec"=>th[1].decimal))
            push!(dict_entry,("bin"=>convert_dec2bin(th[1].decimal)))
            push!(dict_entry,("cost"=>th[1].cost))
            push!(dict_entry,("kV"=>[th[1].PCCcbls[1].elec.volt]))
            push!(dict_entry,("dec_set"=>[th[1].decimal]))
            xfo_cost=0
            conv_cost=0
            for xfo in th[1].pcc.xfmrs
                xfo_cost=xfo_cost+xfo.costs.ttl
            end
            for conv in th[1].pcc.conv
                conv_cost=conv_cost+conv.costs.ttl
            end
            push!(dict_entry,("xfo_cost"=>xfo_cost))
            push!(dict_entry,("conv_cost"=>conv_cost))
            push!(dict_entry,("id"=>th[1].id))
        else
            push!(dict_entry,("mva"=>Float64[]))
            push!(dict_entry,("dec"=>counting_number))
            push!(dict_entry,("bin"=>convert_dec2bin(max_dec)))
            push!(dict_entry,("cost"=>Inf))
            push!(dict_entry,("xfo_cost"=>0.0))
            push!(dict_entry,("conv_cost"=>0.0))
            push!(dict_entry,("kV"=>Float64[]))
            push!(dict_entry,("dec_set"=>Float64[]))
            push!(dict_entry,("id"=>""))
        end
        counting_number=counting_number+1
        push!(dict_of_connections,deepcopy(dict_entry))
    end
    if (max_dec>length(dict_of_connections))
        for i=length(dict_of_connections)+1:1:max_dec
            dict_entry=Dict{String,Any}()
            push!(dict_entry,("mva"=>Float64[]))
            push!(dict_entry,("dec"=>i))
            push!(dict_entry,("bin"=>convert_dec2bin(max_dec)))
            push!(dict_entry,("cost"=>Inf))
            push!(dict_entry,("xfo_cost"=>0.0))
            push!(dict_entry,("conv_cost"=>0.0))
            push!(dict_entry,("kV"=>Float64[]))
            push!(dict_entry,("dec_set"=>Float64[]))
            push!(dict_entry,("id"=>""))
            push!(dict_of_connections,deepcopy(dict_entry))
        end
    end
    return dict_of_connections
end

function combine_TE_parts(child_circ,mom_circ,dad_circ,ocn,ks)
    wnd_power=create_wind_and_power_dict()
    child_circ.oss=[]
    child_circ.mog[1].xfmrs=[]
    child_circ.mog[1].conv=[]
    child_circ.mog[1].plat=[]
    child_circ.cost=0
    child_circ.MVcbls=[]
    child_circ.HVcbls=[]
    child_circ.O2Ocbls=[]
    child_circ.PCCcbls=[]
    child_circ.id=mom_circ.id*"_"*dad_circ.id
    child_circ, wnd_power = transfer_2_circuits_equipment(child_circ,mom_circ,wnd_power)
    child_circ, wnd_power = transfer_2_circuits_equipment(child_circ,dad_circ,wnd_power)
    if (child_circ.mog[1].kV==300.0)
        child_circ.mog[1].conv=set_oss_convertersO2O(child_circ.mog[1].kV,wnd_power,ocn.database,ks)
    end
    #find transformers
    child_circ.mog[1].xfmrs=set_oss_transformers(child_circ.mog[1].kV,wnd_power,ocn.database,ks)
    #set platform structure
    child_circ=set_mog_platformO2O(child_circ,ks)
    child_circ=total_circuit_cost(child_circ)
    return child_circ
end

#copies equipment from circ into child_circ
function transfer_2_circuits_equipment(child_circ,circ,wnd_power)
    for mv_c in circ.MVcbls
        if (mv_c.path[2].num==circ.mog[1].node.num)
            wnd_power=track_wind_and_power(wnd_power,mv_c)
        end
        push!(child_circ.MVcbls,mv_c)
    end
    for hv_c in circ.HVcbls
        if (hv_c.path[2].num==circ.mog[1].node.num)
            wnd_power=track_wind_and_power(wnd_power,hv_c)
        end
        push!(child_circ.HVcbls,hv_c)
    end
    for hv_c in circ.O2Ocbls
        if (hv_c.path[2].num==circ.mog[1].node.num)
            wnd_power=track_wind_and_power(wnd_power,hv_c)
        end
        push!(child_circ.O2Ocbls,hv_c)
    end
    for os in circ.oss
        push!(child_circ.oss,os)
    end
    return child_circ, wnd_power
end

function combine_TH_parts(child_circ,circ,bin,dec,ocn,ks)
    child_circ.binary=bin
    child_circ.decimal=dec
    wnd_power=create_wind_and_power_dict()
    for owp in circ.owpps
        push!(child_circ.owpps,owp)
    end
    for mv_c in circ.MVcbls
        push!(child_circ.MVcbls,mv_c)
    end
    for hv_c in circ.HVcbls
        push!(child_circ.HVcbls,hv_c)
    end
    for hv_c in circ.O2Ocbls
        push!(child_circ.O2Ocbls,hv_c)
    end
    for hv_c in circ.PCCcbls
        push!(child_circ.PCCcbls,hv_c)
    end
    for os in circ.oss
        push!(child_circ.oss,os)
    end
    for os in circ.mog
        push!(child_circ.mog,os)
    end
    for hv_c in child_circ.PCCcbls
        wnd_power=track_wind_and_power(wnd_power,hv_c)
    end
    #NOTE not finished!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    child_circ.pcc.xfmrs=[]
    child_circ.pcc.conv=[]
    if (wnd_power["220.0"]["mva"]!=0.0 && child_circ.pcc.kV!=220.0)
        xfo=deepcopy(ocn.database["transformers"][string(wnd_power["220.0"]["mva"])]["onshore"])
        xfo.wnd=find_netWind(wnd_power["220.0"]["wind"])
        xfo.elec.hv=circ.pcc.kV
        xfo.elec.lv=220.0
        xfo=cost_xfo_pcc(xfo,ks)
        push!(child_circ.pcc.xfmrs,deepcopy(xfo))
    end
    if (wnd_power["400.0"]["mva"]!=0.0 && circ.pcc.kV!=400.0)
        xfo=deepcopy(ocn.database["transformers"][string(wnd_power["400.0"]["mva"])]["onshore"])
        xfo.wnd=find_netWind(wnd_power["400.0"]["wind"])
        xfo.elec.hv=circ.pcc.kV
        xfo.elec.lv=400.0
        xfo=cost_xfo_pcc(xfo,ks)
        push!(child_circ.pcc.xfmrs,deepcopy(xfo))
    end
    if (wnd_power["300.0"]["mva"]!=0.0)
        #converters
        conv=deepcopy(database["converters"][string(wnd_power["300.0"]["mva"])]["onshore"])
        conv.wnd=find_netWind(wnd_power["300.0"]["wind"])
        conv=cost_hvdc_pcc(conv,ks)
        conv=adjust_base_hvdc_onshore_converter(conv,ks)
        push!(child_circ.pcc.conv,deepcopy(conv))

        #transformers
        xfo=deepcopy(database["transformers"][string(wnd_power["300.0"]["mva"])]["onshore"])
        xfo.wnd=find_netWind(wnd_power["300.0"]["wind"])
        xfo.elec.hv=circ.pcc.kV
        xfo.elec.lv=300.0
        xfo=cost_xfo_pcc(xfo,ks)
        push!(child_circ.pcc.xfmrs,deepcopy(xfo))
    end

    child_circ.id=child_circ.id*"_"*circ.id
    child_circ=total_circuit_cost(child_circ)
    return child_circ
end

#te_part =TE_parts[2][1]
function change_tEs_bases_pcc2oss(TE_parts, tc, ocn, ks)
    #changes the connection from the PCC to the MOG
    #the order needs to change find those that total
    for te_parts in TE_parts
        for te_part in te_parts
            te_part=change_tE_base_pcc2oss(deepcopy(tc),te_part,ocn,ks)
        end
    end
    return TE_parts
end

#=
function change_tE_base_pcc2oss(tc,te_part,ocn,ks)
    #foCost=deepcopy(forward.cost)
    km=euclidian_distance(te_part.PCCcbls[1].path[1].xy,tc.mog[1].node.xy)

    if (length(te_part.owpps)>1)
        #find bus to connect
        con_bus=bus()
        wnd_power=create_wind_and_power_dict()
        for os in te_part.mog
            if (os.node.num==te_part.PCCcbls[1].path[1].num)
                con_bus=os
                break
            end
        end

        for mv_c in te_part.MVcbls
            if (mv_c.path[length(mv_c.path)].num==con_bus.node.num)
                wnd_power=track_wind_and_power(wnd_power,mv_c)
            end
        end
        for hv_c in te_part.HVcbls
            if (hv_c.path[length(hv_c.path)].num==con_bus.node.num)
                wnd_power=track_wind_and_power(wnd_power,hv_c)
            end
        end
        for hv_c in te_part.O2Ocbls
            if (hv_c.path[length(hv_c.path)].num==con_bus.node.num)
                wnd_power=track_wind_and_power(wnd_power,hv_c)
            end
        end
        total_power=wnd_power["66.0"]["mva"]+wnd_power["220.0"]["mva"]+wnd_power["400.0"]["mva"]+wnd_power["300.0"]["mva"]
        o2o_options=[oss_acdc_to_oss(km,con_bus,wnd_power,ocn.database,ks)]
        oss_up=te_part.mog[1]
    else#owpp to oss connection
        #Calculate length and type of cables
        te_part.MVcbls=[]
        te_part.HVcbls=[]
        te_part.oss=[]
        total_power=te_part.base.mva
        o2o_options=[owpp_acdc_to_oss(km,te_part.base,tc.mog[1].node,ocn.database,ks,false)]
    end
    te_part.mog=tc.mog
    te_part.PCCcbls=[]
    shore_connection=mog2pcc_possibilities_noMOG(tc.mog[1].node,tc,ocn)
    #find bus to connect
    con_bus2=bus()
    wnd_power=create_wind_and_power_dict()
    o2o_option_dumb=deepcopy(o2o_options[1])
    o2o_option_dumb["300.0"]["cost"]=Inf
    o2o_option_dumb["220.0"]["cost"]=Inf
    o2o_option_dumb["400.0"]["cost"]=Inf
    if (haskey(o2o_option_dumb,"66.0"))
        o2o_option_dumb["66.0"]["cost"]=Inf
    else
        push!(o2o_option_dumb,("66.0"=>deepcopy(o2o_option_dumb["400.0"])))
    end
    for os in tc.mog
        if (os.node.num==tc.PCCcbls[1].path[1].num)
            con_bus2=os
            break
        end
    end

    for mv_c in tc.MVcbls
        if (mv_c.path[length(mv_c.path)].num==con_bus2.node.num)
            o2o_option=deepcopy(o2o_option_dumb)
            o2o_option[string(mv_c.elec.volt)]["cable"]=mv_c
            o2o_option[string(mv_c.elec.volt)]["cost"]=mv_c.costs.grand_ttl
            push!(o2o_options,o2o_option)
        end
    end
    for hv_c in tc.HVcbls
        if (hv_c.path[length(hv_c.path)].num==con_bus2.node.num)
            o2o_option=deepcopy(o2o_option_dumb)
            o2o_option[string(hv_c.elec.volt)]["cable"]=hv_c
            o2o_option[string(hv_c.elec.volt)]["cost"]=hv_c.costs.grand_ttl
            push!(o2o_options,o2o_option)
        end
    end
    for hv_c in tc.O2Ocbls
        if (hv_c.path[length(hv_c.path)].num==con_bus2.node.num)
            o2o_option=deepcopy(o2o_option_dumb)
            o2o_option[string(hv_c.elec.volt)]["cable"]=hv_c
            o2o_option[string(hv_c.elec.volt)]["cost"]=hv_c.costs.grand_ttl
            push!(o2o_options,o2o_option)
        end
    end

    if (length(te_part.owpps)>1)
        te_part=find_optimal_circuitO2O([(km,oss_up)],o2o_options,shore_connection,ocn,te_part,ks)
    else
        println(o2o_options[1]["220.0"]["cable"].mva)
        println(o2o_options[2]["66.0"]["cable"].mva)
        println(shore_connection["220.0"]["cable"].mva)
        te_part=find_optimal_circuit([(km,te_part.base)],o2o_options,shore_connection,ocn,te_part,ks)

        if (tc.PCCcbls[1].elec.volt==300.0)
            if ((length(te_part.HVcbls)==0) || (te_part.HVcbls[1].elec.volt!=300.0))
                te_part.mog[1].plat[2].costs.ttl=te_part.mog[1].plat[2].costs.ttl-(ks.pdc_h+ks.pdc_h*ks.opx_pl*npv_years())
                te_part.mog[1].conv[1].costs.ttl=te_part.mog[1].conv[1].costs.ttl-(ks.conv_d+ks.conv_d*ks.opx_co*npv_years())
            else
                te_part.mog[1].plat[1].costs.ttl=te_part.mog[1].plat[1].costs.ttl-(ks.pac_f+ks.pac_f*ks.opx_pl*npv_years())
            end
        else
            if ((length(te_part.HVcbls)==0) || (te_part.HVcbls[1].elec.volt!=300.0))
                te_part.mog[1].plat[1].costs.ttl=te_part.mog[1].plat[1].costs.ttl-(ks.pac_f+ks.pac_f*ks.opx_pl*npv_years())
            else
            end
        end
    end
    te_part.PCCcbls=[]
    te_part.pcc.xfmrs=[]
    te_part.pcc.conv=[]
    te_part=total_circuit_cost(te_part)
    return te_part
end=#

function find_set_TE_parts(above_beta, circuits_set, above_beta_1s, tc, owpps)
    above_beta_alpha=round(Int32,2^(above_beta_1s[1]-1))#find the alpha+1 bit
    above_beta_0s=findall(x->x==0,above_beta)#find all 0s
    above_beta_dec=round(Int32,convert_bin2dec(above_beta))#find decimal equivalent of binary
    above_beta_dec_no_alpha=above_beta_dec-above_beta_alpha#remove alpha+1 bit
    above_beta_mid_point=round((above_beta_dec/2))#find halfway point for above_beta
    E_parts=make_set_B(owpps,tc.pcc)#B will be reduced to E
    TE_parts=Array{Array{circuit, 1}, 1}()

    mid_point_indx=0
    set=false
    for (i_E, row_E) in enumerate(E_parts)#filter set
        row_E_1s=findall(x->x==1,row_E[1])
        row_E_dec=row_E[2]
        keep=true
        #filter for alpha entry and alpha to gamma inclusive entry
        if (row_E_dec==above_beta_alpha || row_E_dec==above_beta_alpha)
            keep=false
        end
        #check bits in entries are partt of tbj
        if (keep==true)
            for posi_1 in row_E_1s
                #don't keep if a bit=1 when it should equal 0
                for posi_0 in above_beta_0s
                    if (posi_1[2]==posi_0)
                        keep=false
                        break
                    end
                end
                if (keep==false)
                    break
                end
            end
        end
        #Keep circuit equivalent of all remining entries
        if (keep==true)
            if (row_E_dec>above_beta_mid_point && set==false)
                mid_point_indx=deepcopy(length(TE_parts)+1)
                set=true
            end
            push!(TE_parts,deepcopy(circuits_set[i_E]))
        end
    end
    return TE_parts, mid_point_indx
end
################################################################################
function change_tE_base_pcc2oss(tc,te_part,ocn,ks)
    #foCost=deepcopy(forward.cost)
    km=euclidian_distance(te_part.PCCcbls[1].path[1].xy,tc.mog[1].node.xy)
    wnd_power=create_wind_and_power_dict()

    if (length(te_part.owpps)>1)
        #find bus to connect
        con_bus=bus()
        wnd_power=create_wind_and_power_dict()
        for os in te_part.mog
            if (os.node.num==te_part.PCCcbls[1].path[1].num)
                con_bus=os
                break
            end
        end
        for mv_c in te_part.MVcbls
            if (mv_c.path[length(mv_c.path)].num==con_bus.node.num)
                wnd_power=track_wind_and_power(wnd_power,mv_c)
            end
        end
        for hv_c in te_part.HVcbls
            if (hv_c.path[length(hv_c.path)].num==con_bus.node.num)
                wnd_power=track_wind_and_power(wnd_power,hv_c)
            end
        end
        for hv_c in te_part.O2Ocbls
            if (hv_c.path[length(hv_c.path)].num==con_bus.node.num)
                wnd_power=track_wind_and_power(wnd_power,hv_c)
            end
        end
        total_power=wnd_power["66.0"]["mva"]+wnd_power["220.0"]["mva"]+wnd_power["400.0"]["mva"]+wnd_power["300.0"]["mva"]
        o2o_options=[oss_acdc_to_oss(km,con_bus,wnd_power,ocn.database,ks)]
        oss_up=te_part.mog[1]
    else#owpp to oss connection
        #Calculate length and type of cables
        te_part.MVcbls=[]
        te_part.HVcbls=[]
        te_part.oss=[]
        total_power=te_part.base.mva
        o2o_options=[owpp_acdc_to_oss(km,te_part.base,tc.mog[1].node,ocn.database,ks,false)]
    end
    te_part.mog=tc.mog
    te_part.PCCcbls=[]
    shore_connection=Dict{String, Dict{String,Any}}()
    noConKV=Dict{String, Any}()
    conKV=Dict{String, Any}()
    push!(conKV,("cost"=>tc.cost))
    push!(conKV,("cable"=>tc.PCCcbls[1]))
    conKV["cable"].mva=total_power
    xfo_for=transformer()
    xfo_for.costs.ttl=0
    push!(conKV,("xfo_for"=>xfo_for))
    push!(noConKV,("cost"=>Inf))
    if (tc.PCCcbls[1].elec.volt==300.0)
        conv=converter()
        push!(conKV,("conv_for"=>conv))
        conv.costs.ttl=0
        push!(shore_connection,("220.0"=>noConKV))
        push!(shore_connection,("400.0"=>noConKV))
        push!(shore_connection,("300.0"=>conKV))
    elseif (tc.PCCcbls[1].elec.volt==400.0)
        push!(shore_connection,("220.0"=>noConKV))
        push!(shore_connection,("400.0"=>conKV))
        push!(shore_connection,("300.0"=>noConKV))
    elseif (tc.PCCcbls[1].elec.volt==220.0)
        push!(shore_connection,("220.0"=>conKV))
        push!(shore_connection,("400.0"=>noConKV))
        push!(shore_connection,("300.0"=>noConKV))
    end
    if (length(te_part.owpps)>1)
        te_part=find_optimal_circuitO2O([(km,oss_up)],o2o_options,shore_connection,ocn,te_part,ks)
    else
        te_part=find_optimal_circuit([(km,te_part.base)],o2o_options,shore_connection,ocn,te_part,ks)
        #OSS node correction
        #=if (length(te_part.HVcbls)>0)
            ocn.num=ocn.num-1
        end=#
        #plarform and converter corrections
        if (tc.PCCcbls[1].elec.volt==300.0)
            if ((length(te_part.HVcbls)==0) || (te_part.HVcbls[1].elec.volt!=300.0))
                te_part.mog[1].plat[2].costs.ttl=te_part.mog[1].plat[2].costs.ttl-(ks.pdc_h+ks.pdc_h*ks.opx_pl*npv_years())
                te_part.mog[1].conv[1].costs.ttl=te_part.mog[1].conv[1].costs.ttl-(ks.conv_d+ks.conv_d*ks.opx_co*npv_years())
            else
                te_part.mog[1].plat[1].costs.ttl=te_part.mog[1].plat[1].costs.ttl-(ks.pac_f+ks.pac_f*ks.opx_pl*npv_years())
            end
        else
            if ((length(te_part.HVcbls)==0) || (te_part.HVcbls[1].elec.volt!=300.0))
                te_part.mog[1].plat[1].costs.ttl=te_part.mog[1].plat[1].costs.ttl-(ks.pac_f+ks.pac_f*ks.opx_pl*npv_years())
            else
            end
        end
    end
    te_part.PCCcbls=[]
    te_part.pcc.xfmrs=[]
    te_part.pcc.conv=[]
    te_part=total_circuit_cost(te_part)
    return te_part
end
