ks=get_Cost_Data()
df = DataFrame(XLSX.readtable("v2.0/economics/input_data/data.xlsx", "wind_data")...)

conv=converter()

conv.relia=get_onshore_conv_failure_data(conv)
conv.wnd=wndF_wndPrf([getproperty(df,Symbol("A6"))])
conv.mva=1000
conv.elec.mva=1000
conv=cost_hvdc_oss(conv,ks)
#conv=cost_hvdc_pcc(conv,ks)
println("CAPEX: "*string(conv.costs.cpx))
println("OPEX: "*string(conv.costs.cm))
println("CAPEX+OPEX: "*string(conv.costs.cpx+conv.costs.cm))
println("TLC: "*string(conv.costs.tlc))
println("EENS: "*string(conv.costs.eens))
println("Total: "*string(conv.costs.ttl))
println()

#HVDC/HVAC cable testing complete see reference costs Catapult Benchmarking HVDC_cable, HVAC_cable
#HV/MVAC cable testing
save("v2.0/diverse_wind.jld2", "wnd_diverse",wnd_diverse)
save("v2.0/wnd_non_diverse.jld2", "wnd_diverse",wnd_non_diverse)
cbl=cable()
cbl_data=get_220kV_cables()
wnd_diverse=ocean.hv_circuits[length(ocean.hv_circuits)][1].osss_mog[1].wnd
wnd_non_diverse=ocean.hv_circuits[1][1].osss_mog[1].wnd


cbl=cable()
cbl.wnd=wndF_wndPrf([getproperty(df,Symbol("Norther"))])
cbl.mva=700#220kv,400mm,1,25.05km
km=150
num=2
#cbl.wnd.delta=0.38
#fillOut_cable_struct_dc(cbl,cbl_data[1],km,num)
#cost_hvdc_cable(cbl,ks)
#delta → 0.23857749
#lf → 0.37666464

test_hvdc_cable(cbl,cbl_data.cbls150kV,km,ks)
println()
test_hvdc_cable(cbl,cbl_data.cbls300kV,km,ks)
println()

#fillOut_cable_struct_dc(cbl,cbl_data.cbls150kV[1],50,2)
#fillOut_cable_struct_ac(cbl,cbl_data.cbls220kV[1],km,1)
#cbl=cost_mvac_cable(cbl,ks)
#test_mvac_cable(cbl,cbl_data.cbls33kV,km,ks)
test_hvac_cable(cbl,cbl_data,km,ks)
println()
test_hvac_cable(cbl,cbl_data.cbls220kV,km,ks)
println()
test_hvac_cable(cbl,cbl_data.cbls400kV,km,ks)
println()



for mva in [2000,1900,1800,1700,1600,1500,1400,1300,1200,1100,1000]
    plat=platform()
    plat.mva=mva
    plat=cost_dc_platform(plat,ks)
    println(plat.costs.ttl)
end
#transformer testing
xfo=transformer()
xfo_data=get_Xfo_Data()
xfo.relia=get_offshore_xfo_failure_data(xfo)
xfo.wnd=wndF_wndPrf([getproperty(df,Symbol("Norther"))])
xfo.mva=1400
xfo.num=2
xfo.elec.mva=700
#cost_xfo_oss(xfo,ks)
cost_xfo_oss(xfo,ks)
println()

plat=platform()
plat.mva=100
cost_ac_platform(plat,ks)
println()
println("MVA: "*string(plat.mva))
println("CAPEX: "*string(plat.costs.cpx))
println("OPEX: "*string(plat.costs.cm))
println("Total: "*string(plat.costs.ttl))

xfo=transformer()
xfo.relia=get_onshore_xfo_failure_data(xfo)
xfo.wnd=wnd
xfo.mva=700
test_xfo_pcc(xfo,ks)
println()
#xfo=cost_xfo_pcc(xfo,ks)

function test_hvdc_cable(cbl0,cds,km,ks)
    cbl=cable()
    cbl.costs.ttl=Inf
    for cd in cds
        num=2
        while (cbl0.mva>(num*cd[1]*cd[5])*10^-3)
            num=num+2
        end
        fillOut_cable_struct_dc(cbl0,cd,km,num)
        cbl0=cost_hvdc_cable(cbl0,ks)
        if (cbl0.costs.ttl<cbl.costs.ttl)
            cbl=deepcopy(cbl0)
        end

    end
    #println("Cable: "*string(cbl.elec.volt)*" - "*string(cbl.num)*" - "*string(cbl.size)*" Cable MVA: "*string(cbl.elec.mva*cbl.num)*" Circuit MVA: "*string(cbl.mva)*", km: "*string(cbl.length))
    println("Cable: "*string(cbl.elec.volt)*" - "*string(cbl.num)*" - "*string(cbl.size))
    println("PROCUREMENT: "*string(cbl.costs.cpx_p))
    println("INSTALLATION: "*string(cbl.costs.cpx_i))
    println("CAPEX: "*string(cbl.costs.cpx_p+cbl.costs.cpx_i))
    println("OPEX: "*string(cbl.costs.cm))
    println("CAPEX+OPEX: "*string(cbl.costs.cpx_p+cbl.costs.cpx_i+cbl.costs.cm))
    println("LOSSES: "*string(cbl.costs.rlc))
    println("EENS: "*string(cbl.costs.eens))
    println("Total: "*string(cbl.costs.ttl))
end

function test_hvac_cable(cbl0,cds,km,ks)
    cbl=cable()
    cbl.costs.ttl=Inf
    for cd in cds
        num=1
        while (cbl0.mva>num*get_newQ_Capacity(cbl0.elec.freq,km,cd[1],cd[4],cd[5]))
            num=deepcopy(num+1)
        end
        fillOut_cable_struct_ac(cbl0,cd,km,num)
        #cbl0=cost_hvac_cable_o2o(cbl0,ks)
        cbl0=cost_hvac_cable(cbl0,ks)
        if (cbl0.costs.ttl<cbl.costs.ttl)
            cbl=deepcopy(cbl0)
        end
    end
    #println("Cable: "*string(cbl.elec.volt)*" - "*string(cbl.num)*" - "*string(cbl.size)*" Cable MVA: "*string(cbl.elec.mva*cbl.num)*" Circuit MVA: "*string(cbl.mva)*", km: "*string(cbl.length))
    println("Cable: "*string(cbl.elec.volt)*" - "*string(cbl.num)*" - "*string(cbl.size))
    println("PROCUREMENT: "*string(cbl.costs.cpx_p))
    println("INSTALLATION: "*string(cbl.costs.cpx_i))
    println("CAPEX: "*string(cbl.costs.cpx_p+cbl.costs.cpx_i))
    println("OPEX: "*string(cbl.costs.cm))
    println("CAPEX+OPEX: "*string(cbl.costs.cpx_p+cbl.costs.cpx_i+cbl.costs.cm))
    println("QC: "*string(cbl.costs.qc))
    println("REACTORS: "*string(cbl.reactors))
    println("S.G.: "*string(cbl.costs.sg))
    println("LOSSES: "*string(cbl.costs.rlc))
    println("EENS: "*string(cbl.costs.eens))
    println("Total: "*string(cbl.costs.ttl))
end

function test_mvac_cable(cbl0,cds,km,ks)
    cbl=cable()
    cbl.costs.ttl=Inf
    for cd in cds
        num=1
        while (cbl0.mva>num*get_newQ_Capacity(cbl0.elec.freq,km,cd[1],cd[4],cd[5]))
            num=deepcopy(num+1)
        end
        fillOut_cable_struct_ac(cbl0,cd,km,num)
        #cbl0=cost_hvac_cable_o2o(cbl0,ks)
        cbl0=cost_mvac_cable(cbl0,ks)
        if (cbl0.costs.ttl<cbl.costs.ttl)
            cbl=deepcopy(cbl0)
        end
    end
    println("Cable: "*string(cbl.elec.volt)*" - "*string(cbl.num)*" - "*string(cbl.size)*" Cable MVA: "*string(cbl.elec.mva*cbl.num)*" Circuit MVA: "*string(cbl.mva)*", km: "*string(cbl.length))
    println("PROCUREMENT: "*string(cbl.costs.cpx_p))
    println("INSTALLATION: "*string(cbl.costs.cpx_i))
    println("CAPEX: "*string(cbl.costs.cpx_p+cbl.costs.cpx_i))
    println("OPEX: "*string(cbl.costs.cm))
    println("CAPEX+OPEX: "*string(cbl.costs.cpx_p+cbl.costs.cpx_i+cbl.costs.cm))
    println("QC: "*string(cbl.costs.qc))
    println("LOSSES: "*string(cbl.costs.rlc))
    println("EENS: "*string(cbl.costs.eens))
    println("Total: "*string(cbl.costs.ttl))
end

function test_xfo_oss(xfo0,ks)
    xfo=transformer()
    xfo.costs.ttl=Inf
    xfo_data=get_Xfo_Data()
    for xd in xfo_data
        if ((xd+10<=xfo0.mva) && (3*xd>xfo0.mva))
            xfo0.num=1
            xfo0.elec.mva=xd
            while ((xfo0.num*xfo0.elec.mva)<xfo0.mva)
                xfo0.num=xfo0.num+1
            end
            xfo0=cost_xfo_oss(xfo0,ks)
            if (xfo0.costs.ttl<xfo.costs.ttl)
                xfo=deepcopy(xfo0)
            end
        end
    end
    println("XFO: "*string(xfo.elec.mva)*" x "*string(xfo.num))
    println("PROCUREMENT: "*string(xfo.costs.cpx_p))
    println("INSTALLATION: "*string(xfo.costs.cpx_i))
    println("CAPEX: "*string(xfo.costs.cpx_p+xfo.costs.cpx_i))
    println("OPEX: "*string(xfo.costs.cm))
    println("CAPEX+OPEX: "*string(xfo.costs.cpx_p+xfo.costs.cpx_i+xfo.costs.cm))
    println("LOSSES: "*string(xfo.costs.tlc))
    println("EENS: "*string(xfo.costs.eens))
    println("Total: "*string(xfo.costs.ttl))
end


function test_xfo_pcc(xfo0,ks)
    xfo=transformer()
    xfo.costs.ttl=Inf
    xfo_data=get_Xfo_Data()
    for xd in xfo_data
        if ((xd+10<=xfo0.mva) && (3*xd>xfo0.mva))
            xfo0.num=1
            xfo0.elec.mva=xd
            while ((xfo0.num*xfo0.elec.mva)<xfo0.mva)
                xfo0.num=xfo0.num+1
            end
            xfo0=cost_xfo_pcc(xfo0,ks)
            if (xfo0.costs.ttl<xfo.costs.ttl)
                xfo=deepcopy(xfo0)
            end
        end
    end
    println("XFO: "*string(xfo.elec.mva)*" x "*string(xfo.num))
    println("PROCUREMENT: "*string(xfo.costs.cpx_p))
    println("INSTALLATION: "*string(xfo.costs.cpx_i))
    println("CAPEX: "*string(xfo.costs.cpx_p+xfo.costs.cpx_i))
    println("OPEX: "*string(xfo.costs.cm))
    println("CAPEX+OPEX: "*string(xfo.costs.cpx_p+xfo.costs.cpx_i+xfo.costs.cm))
    println("LOSSES: "*string(xfo.costs.tlc))
    println("EENS: "*string(xfo.costs.eens))
    println("Total: "*string(xfo.costs.ttl))
end


#################### Depricated

#=
#calculates the cost of a given hvac cable OSS to PCC
function cost_hvac_cable_o2p(cbl,ks)
    #cost of losses in the cable
    cbl.costs.rlc=cost_rlc(cbl,ks)
    #cost of cable compensation placed on OSS - divide by 2 for each OSS
    cbl.costs.qc_oss,cbl.costs.qc_pcc=cost_qc_hvac_pcc(cbl,ks)
    #capex of cable
    cbl.costs.cpx_p,cbl.costs.cpx_i=capex_hvac_cable(cbl,ks)
    #cost of corrective maintenance
    cbl.costs.cm=cost_cm(cbl,ks)
    #cost of expected energy not served
    cbl.costs.eens=cost_eens(cbl,ks)
    #totals the cable cost
    cbl.costs.ttl,cbl.costs.perkm_ttl=cost_cbl_sum(cbl)
    return cbl
end
=#


#=
#HVAC compensation cost oss 2 pcc
function cost_qc_hvac_pcc(cbl,ks)
    #div sets compensation to 50-50 split
    div=0.5
    f=cbl.elec.freq
    A=cbl.elec.farrad*cbl.length*cbl.num
    Q=2*pi*f*cbl.elec.volt^2*A
    Q_oss=Q*div
    Q_pcc=Q*(1-div)
    qc_oss=ks.Qc_oss*Q_oss
    qc_pcc=ks.Qc_pcc*Q_pcc
    return qc_oss, qc_pcc
end=#
