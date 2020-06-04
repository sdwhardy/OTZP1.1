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

function copy_below_beta_equipment(crc,beta_n_below,ocn,beta)
    removedMVmain=false
    oss_system=circuit()
    #base_bin=top_dec2bin(base_dec)
    #base quantities
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

#oss=1526

function add_temporary_hv_connections(forwardMV,ocn,beta)
    ks=get_Cost_Data()
    destination=forwardMV.mog[1].node
    #=for mv_c in forwardMV.MVcbls
        if (mv_c.path[length(mv_c.path)].num==destination.num)
             mv_c.path[length(mv_c.path)]=destination
        end
    end=#

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
