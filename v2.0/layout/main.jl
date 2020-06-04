include("functions.jl")

#pcc,owpps,offset=utm_gps2xy(pcs[2],owps,32,true)
#main ocean layout function
function ocean_layout(pcc,owpps,offset)
    println("Laying out the Offshore Region...")
    ocean=eez()#builds eez
    ocean.offset=offset#saves the offset
    ocean.pcc=pcc
    ocean.owpps=owpps
    ocean.owpps=order_owpps(ocean.owpps,pcc)
    ocean=number_buses(ocean)
    ocean.sys=get_SystemData()
    println("Building Equipment Reference Table...")
    km_mva_set=length_power_set(ocean.owpps)
    ocean.database=get_equipment_tables(km_mva_set,ocean.owpps[1].wnd,pcc.kV)
    println("Finding Max MV Range...")
    ocean.owpps=find_max_mv_transmission(ocean.owpps,ocean.database)
    println("Finding Basis Set Tb...")
    @time B=make_set_B(ocean.owpps,pcc)
    @time ocean=make_set_Tb(ocean.owpps,ocean.pcc,B,ocean)
    println("Optimizing Tb CIrcuits...")
    save("v2.0/tempfiles/ocean/b4Opt_tnep_ocean.jld2", "ocean",ocean)
    @time ocean.hv_circuits=finalize_circuits_layout(ocean.hv_circuits,ocean)
    @time ocean.mv_circuits=finalize_circuits_layout(ocean.mv_circuits,ocean)
    save("v2.0/tempfiles/ocean/greedy_tnep_ocean.jld2", "ocean",ocean)
    return ocean
end

#finds the maximum length required for any single power level
function length_power_set(owpps)
    power_pairs=Array{Tuple,1}()
    for wp in owpps
        push!(power_pairs,(wp.node.num,wp.mva))
    end
    power_pair_sets=collect(combinations(power_pairs))
    power_pairs=Array{Tuple,1}()
    for power_pair_set in power_pair_sets
        base_wp=Inf
        power_sum=0
        for power_pair in power_pair_set
            power_sum=power_sum+power_pair[2]
            if (base_wp>power_pair[1])
                base_wp=deepcopy(power_pair[1])
            end
        end
        push!(power_pairs,(base_wp,power_sum))
    end
    unique!(power_pairs)
    sort!(power_pairs, by = x -> x[2],rev=false)
    km_mva_set=Array{Tuple,1}()
    push!(km_mva_set,power_pairs[1])
    for power_pair in power_pairs[2:length(power_pairs)]
        if (km_mva_set[length(km_mva_set)][2]==power_pair[2] && km_mva_set[length(km_mva_set)][1]<power_pair[1])
            km_mva_set[length(km_mva_set)]=deepcopy(power_pair)
        elseif (km_mva_set[length(km_mva_set)][2]<power_pair[2])
            push!(km_mva_set,deepcopy(power_pair))
        end
    end
    #=for (index,length_power) in enumerate(km_mva_set)
        km=owpps[length_power[1]].mv_zone
        km_mva_set[index]=(km*1.1,length_power[2])#gives 10% margin of upper limit for cable sizes in look up table
    end=#
    for (index,length_power) in enumerate(km_mva_set)
        km=owpps[length(owpps)].mv_zone
        km_mva_set[index]=(km,length_power[2])#sets longest range for each power level
    end
    return km_mva_set
end
