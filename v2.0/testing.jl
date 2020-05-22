pcs=get_PccData()
owps=get_OwppData()
pccs=pcs
owpps=owps
pcc,owpps,offset=utm_gps2xy(pccs[2],owpps,32,true)
ocean=eez()#builds eez
ocean.offset=offset#saves the offset
ocean.pcc=pcc
ocean.owpps=owpps
ocean.owpps=order_owpps(ocean.owpps,pcc)
ocean=number_buses(ocean)
ocean.sys=get_SystemData()
km_mva_set=length_power_set(ocean.owpps)
ocean.database=get_equipment_tables(km_mva_set,ocean.owpps[1].wnd,220.0)
ocean.owpps=find_max_mv_transmission(ocean.owpps,ocean.database)
con=optimal_owpp2pcc(ocean.owpps[1],pcc,ocean.database)
km=euclidian_distance(ocean.owpps[1].node.xy,pcc.node.xy)
ocean.owpps[1].mva=2500.0
ocean.owpps[1].node.xy.y=240
km_mva_set[1]=(250.0,7500.0)
km_mva_set[2]=(250.0,10000.0)
km_mva_set[3]=(250.0,1900.0)
km_mva_set[4]=(250.0,2000.0)
km_mva_set[5]=(250.0,2500.0)
km_mva_set[6]=(250.0,3000.0)
km_mva_set[7]=(250.0,4000.0)
km_mva_set[8]=(250.0,5000.0)
ks=get_Cost_Data()
conv=ocean.database["converters"]["600.0"]["offshore"]
conv=adjust_base_hvdc_offshore_converter(conv,ks)

#150km and 1400MW should match result in reference costs excel HV test sheet
dc1400=hvdc_canditate_mog2pcc(190.0,1500.0,220.0,ocean.owpps[1].wnd,ocean.database,ks)
print_hvdc_canditate_mog2pcc_result(dc1400)
println()
ac1400=hvac_canditate_mog2pcc(190.0,1500.0,220.0,ocean.owpps[1].wnd,ocean.database,ks)
print_hvac_canditate_mog2pcc_result(ac1400)
#println("CAPEX ttl: "dc1200["cable"].costs.cpx_p+dc1200["cable"].costs.cpx_i+2*dc1200["conv_for"].costs.cpx+2*dc1200["xfo_for"].costs.cpx_p+2*dc1200["xfo_for"].costs.cpx_i)

optimal_mog2pcc(100.0,650.0,220.0,ocean.owpps[1].wnd,ocean.database)
for mv in ["250.0","500.0","750.0","1000.0","1250.0","1500.0","1750.0","2000.0"]
    for c in ocean.database["cables"]["66.0"][mv]
        println(c.mva)
        println(c.size)
        println(c.num)
    end
end
for l=50:10:250
    c=optimal_mog2pcc(l,10000.0,220.0,ocean.owpps[1].wnd,ocean.database)
    if (c["cable"].mpc_ac==true)
        print("MPC - ")
    elseif (c["cable"].elec.volt==300.0)
        print("HVDC - ")
    else
        print("HVAC - ")
    end
    println("TTL: "*string(c["cost"]))
end
println()


#NOTE Below still untested
owpp2pcc=optimal_owpp2pcc(ocn.owpps[1],12,ocn.pcc.kV,ocn.database)
owpp2oss=optimal_owpp2oss(ocn.owpps[1],12,ocn.pcc.kV,ocn.database)

#Testing below
plot_PCCnOWPP(ocn)

println(km_mva_set)

km_mva_set=length_power_set(ocn.owpps)
#for kv in ["220.0","400.0","66.0","150.0","300.0"]
for kv in ["220.0","400.0","66.0","300.0"]
    for mva in km_mva_set
        for cbl in ocn.database["cables"][kv][string(mva[2])]
            println("circuit: "*string(mva[2])*"cable: "*string(cbl.mva)*" - "*string(cbl.elec.volt)*" - "*string(cbl.num)*" - "*string(cbl.size)*" @ "*string(cbl.length)*" ME "*string(cbl.costs.ttl))
        end
    end
end
println()
for mva in km_mva_set
    for location in ["offshore","onshore"]
        println("circuit: "*string(mva[2])*"xfo: "*string(ocn.database["transformers"][string(mva[2])][location].elec.mva)*" - "*" - "*string(ocn.database["transformers"][string(mva[2])][location].num)*" @ "*string(ocn.database["transformers"][string(mva[2])][location].costs.ttl))
    end
end
println()
for mva in km_mva_set
    for location in ["offshore","onshore"]
        println("circuit: "*string(mva[2])*"conv: "*string(ocn.database["converters"][string(mva[2])][location].elec.mva)*" - "*" - "*string(ocn.database["converters"][string(mva[2])][location].num)*" @ "*string(ocn.database["converters"][string(mva[2])][location].costs.ttl))
    end
end
