
############################ Reads Data from Excel ########################
################################################################################
function get_Cost_Data()
    ks=cost_ks()
    ks.FC_ac=5.8#fixed AC cost
    ks.FC_dc=29#fixed DC cost
    ks.dc=0.25#penalization factor for different than 2 xfrms
    ks.f_ct=0.0085#generating plant variable cost
    ks.p_ct=0.0085#substructure variable cost
    ks.c_ct=0.1276#hvdc converter variable cost
    ks.Qc_oss=0.029#M£/MVAr
    ks.Qc_pcc=0.0174#M£/MVAr
    ks.life=25#lifetime of wind farm
    ks.T_op=219150#Operational lifetime in hours
    ks.E_op=0.00005#Energy price £/Wh
    ks.cf=30#Capitalization factor
    ks.FC_bld=10#Build cost
    ks.p2e=1.16#pounds/euro exchange
    #the cost of laying cables varies greatly - the original number is directly from North grid, but entso-e's number is less than half this, DTU tech catalogue is around whats given
    #the value for HVAC is half of highest end in North Sea, MVAC is half of lowest end - the dat ais very uncertain.
    #HVDC is per cable therefore is *2, while AC is for a 3 phase cable
    ks.ldc=1.775/2#Cost of laying a single core DC cable
    ks.lac=2.12/2#Cost of laying a single 3 phase AC cable
    ks.lmac=1.69/2#Cost of laying a single 3 phase AC cable
    #ks.lmac=1.905/2#Cost of laying a single 3 phase AC cable
    ks.npv=1.04#discount rate - not used calc is done inline change at npv_hours(), npv_years()
    ks.opx_c=0.025#OPEX percent
    ks.opx_pl=0.02#dc cable
    ks.opx_x=0.0015#Offshore transformer
    ks.opx_co=0.02#onshore transformer
    ks.opx_cp=0.007#ac platform
    ks.conv_c=0.0589#c value of converter cost
    ks.conv_d=54.895#fixed construction converter cost
    ks.pac_e=0.105#e value of ac platform cost
    ks.pac_f=16#f value of ac platform cost
    ks.pdc_g=0.125#g value of dc platform cost
    ks.pdc_h=165#h value of dc platform cost
    return ks
end
###########################################################
######################### cables ##########################
###########################################################
#original data is in excel file but reading file and can access directly or to put in a struct togther use with get_Cable_Data_from_file()
function get_66kV_cables()
        cbls66=[[66.0, 95.0, 218.0, 0.17, 300.0, 187.92, 0.44], [66.0, 120.0, 172.0, 0.18, 340.0, 199.52, 0.43], [66.0, 150.0, 136.0, 0.19, 375.0, 211.12, 0.41], [66.0, 185.0,
    110.0, 0.2, 420.0, 227.36, 0.4], [66.0, 240.0, 84.8, 0.22, 480.0, 251.72, 0.38], [66.0, 300.0, 67.6, 0.24, 530.0, 274.92, 0.37], [66.0, 400.0, 53.2, 0.26, 590.0, 306.24, 0.35],
    [66.0, 500.0, 42.8, 0.29, 655.0, 345.68, 0.34], [66.0, 630.0, 34.6, 0.32, 715.0, 387.44, 0.33], [66.0, 800.0, 28.7, 0.35, 775.0, 436.16, 0.32], [66.0, 1000.0, 24.5, 0.38, 825.0, 482.56, 0.31]]
    return cbls66
end

function get_220kV_cables()
    cbls220=[[220.0, 400.0, 60.1, 122.0, 665.0, 496.48, 0.457], [220.0, 500.0, 48.9, 136.0, 732.0, 597.4, 0.437], [220.0, 630.0, 39.1, 151.0, 808.0, 638.0, 0.415], [220.0, 800.0, 31.9, 163.0, 879.0, 783.0, 0.4], [220.0, 1000.0, 27.0, 177.0, 942.0, 812.0, 0.386]]
    return cbls220
end

function get_400kV_cables()
        cbls400=[[400.0, 500.0, 40.0, 117.0, 776.0, 828.24, 0.589], [400.0, 630.0, 36.0, 125.0, 824.0, 916.4, 0.561], [400.0, 800.0, 31.4, 130.0, 870.0, 997.6, 0.54], [400.0, 1000.0, 26.5, 140.0, 932.0, 1154.2, 0.52], [400.0, 1200.0, 22.1, 170.0, 986.0, 1310.8, 0.49], [400.0, 1400.0, 18.9, 180.0, 1015.0, 1467.4, 0.47], [400.0, 1600.0, 16.6, 190.0, 1036.0, 1624.0, 0.46], [400.0, 2000.0, 13.2, 200.0, 1078.0, 1780.6, 0.44]]
        #cbls400=[[400.0, 2000.0, 13.2, 200.0, Inf, Inf, Inf]]
        return cbls400
end

function get_150kV_cables()
    cbls150=[[150.0, 1000.0, 22.4, 0.7, 1644.0, 179.8, 0.7], [150.0, 1200.0, 19.2, 0.7, 1791.0, 208.8, 0.7], [150.0, 1400.0, 16.5, 0.7, 1962.0, 234.9, 0.7], [150.0, 1600.0, 14.4, 0.7, 2123.0, 261.0, 0.7], [150.0, 2000.0, 11.5, 0.7, 2407.0, 290.0, 0.7]]
        return cbls150
end

function get_300kV_cables()
    cbls300=[[300.0, 1000.0, 22.4, 0.7, 1644.0, 263.9, 0.7], [300.0, 1200.0, 19.2, 0.7, 1791.0, 295.8, 0.7], [300.0, 1400.0, 16.5, 0.7, 1962.0, 333.5, 0.7], [300.0, 1600.0, 14.4, 0.7, 2123.0, 371.2, 0.7], [300.0, 2000.0, 11.5, 0.7, 2407.0, 411.8, 0.7]]
    return cbls300
end

#Fills in the physical data of a cable into the cable structure
#common data between ac and dc cables - note in dc cables capacitive and inductive reactance is not appliccable
function fillOut_cable_struct(cbl,cbl_data,km,num)
    cbl.elec.volt=cbl_data[1]
    cbl.size=cbl_data[2]
    cbl.elec.ohm=cbl_data[3]*10^-3
    cbl.elec.farrad=cbl_data[4]*10^-9
    cbl.elec.amp=cbl_data[5]
    cbl.costs.perkm_cpx=cbl_data[6]*10^-3
    cbl.length=km
    cbl.elec.henry=cbl_data[7]*10^-3
    cbl.elec.xl=2*pi*cbl.elec.freq*cbl.elec.henry#cable inductive reactance
    cbl.elec.yc=2*pi*cbl.elec.freq*cbl.elec.farrad#cable capacitive reactance
    cbl.num=num
    cbl=get_cbl_failure_data(cbl)#Set failure data
    return cbl
end

#Fills in the physical data of an ac cable into the cable structure **
function fillOut_cable_struct_ac(cbl,cbl_data,km,num)
    cbl=fillOut_cable_struct(cbl,cbl_data,km,num)#common data between ac and dc cables
    cbl.elec.mva=get_newQ_Capacity(cbl.elec.freq,km,cbl.elec.volt,cbl.elec.farrad,cbl.elec.amp)
    return cbl
end

#Fills in the physical data of an hvdc cable into the cable structure **
function fillOut_cable_struct_dc(cbl,cbl_data,km,num)
    cbl=fillOut_cable_struct(cbl,cbl_data,km,num)#common data between ac and dc cables
    cbl.elec.mva=cbl.elec.volt*cbl.elec.amp*10^-3#dc does not change with distance
    return cbl
end


function get_cbl_failure_data(cbl)
    cbl.relia.fr=(0.1114/100)*cbl.length#/yr/100km
    cbl.relia.mttr=1.971#months
    cbl.relia.mc=0.56#No longer used
    return cbl
end

#Calculates the new hvac cable capacity after 50-50 compensation at distance km. **
function get_newQ_Capacity(f,l,v,q,a)
#Calculates the square of new hvac cable capacity after 50-50 compensation at distance km.
    mva=(sqrt(3)*v*10^3*a/10^6)^2-((0.5*((v*10^3)^2*2*pi*f*l*q))/10^6)^2
#takes square root if negative returns zero if negative
    if mva>=0
        mva=sqrt(mva)
    else
        mva=0.0
    end
 return mva
end

###########################################################
################### transformers ##########################
###########################################################
#Sets all options for transformer sizes in 10MVA steps **
#Max size considered 5000MVA => although not realistic it's to test 10GW HVDC
function  get_Xfo_Data()
    xfos=Array{Float32,1}()
    push!(xfos,10)
    for i=50:10:1500
        push!(xfos,i)
    end
    push!(xfos,2500)
    push!(xfos,5000)
    return xfos
end

#failure data for transformers **
function get_offshore_xfo_failure_data(xfo)
    #failure rate is faily consistent among sources but MTTR varies alot
    #Reference from North Sea Grid Annex - 120hrs
    #6 months is most commnly assumed value but seems old, more recently it's dropped to 2 or 3
    #Energy Transmission and Grid Integration of AC Offshore Wind Farms
    #Integration of Large Scale Wind Energy with Electrical Power Systems in China
    #Reliability Study Analysis of Electrical Systems within Offshore Wind Parks Elforsk report 07:65
    xfo.relia.fr=0.02#/yr# from CIGRE
    xfo.relia.mttr=3#month - assumed as it's middle of the values
    #xfo.relia.mc=2.8#
    #xfo.relia.mc=0.05#Not directly sourced but derived from Operation and maintenance requirements (No longer used)
    return xfo.relia
end

#failure data for transformers **
function get_onshore_xfo_failure_data(xfo)
    #failure rate is faily consistent among sources but MTTR varies alot
    #Reference from North Sea Grid Annex - 120hrs
    #6 months is most commnly assumed value but seems old, more recently it's dropped to 2 or 3
    #Energy Transmission and Grid Integration of AC Offshore Wind Farms
    #Integration of Large Scale Wind Energy with Electrical Power Systems in China
    #Reliability Study Analysis of Electrical Systems within Offshore Wind Parks Elforsk report 07:65
    xfo.relia.fr=0.02#/yr# from CIGRE
    xfo.relia.mttr=0.6#month - assumed as it's middle of the values
    #xfo.relia.mc=2.8#
    #xfo.relia.mc=0.05#Not directly sourced but derived from Operation and maintenance requirements (No longer used)
    return xfo.relia
end

#for now both the same for onshore as offshore until etter data is obtained


function get_onshore_conv_failure_data(conv)
    conv.relia.fr=0.12#/yr
    conv.relia.mttr=2#month
    conv.relia.mc=0.56#
    return conv.relia
end

function get_offshore_conv_failure_data(conv)
    conv.relia.fr=0.12#/yr
    conv.relia.mttr=2#month
    conv.relia.mc=0.56#
    return conv.relia
end


#sizes and costs reactors - same reactors each side
#set of reactors
#All reactor sizes and cost from North Sea Grid Final report AnnexAnnex
function get_reactors(Q)
    reactors=Array{Int32,1}()
    qc=0
    while (Q>275)
        push!(reactors,350)
        qc=qc+10.065
        Q=Q-350
    end
    while (Q>150)
        push!(reactors,200)
        qc=qc+5.75
        Q=Q-200
    end
    while (Q>125)
        push!(reactors,150)
        qc=qc+4.315
        Q=Q-150
    end
    while (Q>75)
        push!(reactors,125)
        qc=qc+3.595
        Q=Q-125
    end
    while (Q>50)
        push!(reactors,75)
        qc=qc+2.155
        Q=Q-75
    end
    while (Q>0)
        push!(reactors,50)
        qc=qc+1.437
        Q=Q-50
    end
    if ((length(reactors)==0) && (Q>0 && Q<50))
        push!(reactors,50)
        qc=qc+1.437
    elseif (Q>0)
        println("Get Reactor Function, Problem assigning reactors.")
    end
    return qc,reactors
end


#Gets the data for equipment to be used, not used unless explicitly called, only for updating to new data if desired
function get_Cable_Data_from_file()
    df = DataFrame(XLSX.readtable("v2.0/economics/input_data/cost_data.xlsx", "eqp_data")...)
    eqps=cbls_data()#create an instance of object ks
    for index=1:length(df[:, 1])
        #get 33kv cable data
        if (typeof(df.kv_33[index]) != Missing)
            cbd=Float32[]
            push!(cbd,df.kv_33[index])
            push!(cbd,df.mm2_33[index])
            push!(cbd,df.ohms_per_km_33[index])
            push!(cbd,df.nf_per_km_33[index])
            push!(cbd,df.Amps_33[index])
            push!(cbd,df.euro_per_m_33[index])
            push!(cbd,df.mh_per_km_33[index])
            push!(eqps.cbls33kV,deepcopy(cbd))
        end

        #get 66kv cable data
        if (typeof(df.kv_66[index]) != Missing)
            cbd=Float32[]
            push!(cbd,df.kv_66[index])
            push!(cbd,df.mm2_66[index])
            push!(cbd,df.ohms_per_km_66[index])
            push!(cbd,df.nf_per_km_66[index])
            push!(cbd,df.Amps_66[index])
            push!(cbd,df.euro_per_m_66[index])
            push!(cbd,df.mh_per_km_66[index])
            push!(eqps.cbls66kV,deepcopy(cbd))
        end

        #get 132kv cable data
        if (typeof(df.kv_132[index]) != Missing)
            cbd=Float32[]
            push!(cbd,df.kv_132[index])
            push!(cbd,df.mm2_132[index])
            push!(cbd,df.ohms_per_km_132[index])
            push!(cbd,df.nf_per_km_132[index])
            push!(cbd,df.Amps_132[index])
            push!(cbd,df.euro_per_m_132[index])
            push!(cbd,df.mh_per_km_132[index])
            push!(eqps.cbls132kV,deepcopy(cbd))
        end

        #get 220kv cable data
        if (typeof(df.kv_220[index]) != Missing)
            cbd=Float32[]
            push!(cbd,df.kv_220[index])
            push!(cbd,df.mm2_220[index])
            push!(cbd,df.ohms_per_km_220[index])
            push!(cbd,df.nf_per_km_220[index])
            push!(cbd,df.Amps_220[index])
            push!(cbd,df.euro_per_m_220[index])
            push!(cbd,df.mh_per_km_220[index])
            push!(eqps.cbls220kV,deepcopy(cbd))
        end


        #get 400kv cable data
        if (typeof(df.kv_400[index]) != Missing)
            cbd=Float32[]
            push!(cbd,df.kv_400[index])
            push!(cbd,df.mm2_400[index])
            push!(cbd,df.ohms_per_km_400[index])
            push!(cbd,df.nf_per_km_400[index])
            push!(cbd,df.Amps_400[index])
            push!(cbd,df.euro_per_m_400[index])
            push!(cbd,df.mh_per_km_400[index])
            push!(eqps.cbls400kV,deepcopy(cbd))
        end

        #get 150kv cable data
        if (typeof(df.kv_150[index]) != Missing)
            cbd=Float32[]
            push!(cbd,df.kv_150[index])
            push!(cbd,df.mm2_150[index])
            push!(cbd,df.ohms_per_km_150[index])
            push!(cbd,df.nf_per_km_150[index])
            push!(cbd,df.Amps_150[index])
            push!(cbd,df.euro_per_m_150[index])
            push!(cbd,df.mh_per_km_150[index])
            push!(eqps.cbls150kV,deepcopy(cbd))
        end


        #get 300kv cable data
        if (typeof(df.kv_300[index]) != Missing)
            cbd=Float32[]
            push!(cbd,df.kv_300[index])
            push!(cbd,df.mm2_300[index])
            push!(cbd,df.ohms_per_km_300[index])
            push!(cbd,df.nf_per_km_300[index])
            push!(cbd,df.Amps_300[index])
            push!(cbd,df.euro_per_m_300[index])
            push!(cbd,df.mh_per_km_300[index])
            push!(eqps.cbls300kV,deepcopy(cbd))
        end
    end
    return eqps
end
################################################################################
############################ Depricated ########################################
#was used to read data from file but now its all inline
#Sets values of all cost factors discribed below
#the structure of object ks is described in file cst_structure.jl
#=function get_Cost_Data_file()
    sht = XLSX.readxlsx("v2.0/economics/input_data/cost_data.xlsx")["cost_data"]
    ks=cost_ks()#create an instance of object ks
    ks.FC_ac=sht["B1"]#fixed AC cost
    ks.FC_dc=sht["B2"]#fixed DC cost
    ks.dc=sht["B3"]#penalization factor for different than 2 xfrms
    ks.f_ct=sht["B4"]#generating plant variable cost
    ks.p_ct=sht["B5"]#substructure variable cost
    ks.c_ct=sht["B6"]#hvdc converter variable cost
    ks.Qc_oss=sht["B7"]#M£/MVAr
    ks.Qc_pcc=sht["B8"]#M£/MVAr
    ks.life=sht["B9"]#lifetime of wind farm
    ks.T_op=sht["B10"]#Operational lifetime in hours
    ks.E_op=sht["B11"]#Energy price £/Wh
    ks.cf=sht["B12"]#Capitalization factor
    ks.FC_bld=sht["B13"]#Build cost
    ks.p2e=sht["B14"]#pounds/euro exchange
    ks.ldc=sht["B15"]#Cost of laying a single core DC cable
    ks.lac=sht["B16"]#Cost of laying a single 3 phase AC cable
    ks.lmac=sht["B17"]#Cost of laying a single 3 phase AC cable
    ks.npv=sht["B18"]#discount rate - not used calc is done inline change at npv_hours(), npv_years()
    ks.opx_c=sht["B19"]#OPEX percent: AC cable
    ks.opx_pl=sht["B20"]#dc cable
    ks.opx_x=sht["B21"]#Offshore transformer
    ks.opx_co=sht["B22"]#onshore transformer
    ks.opx_cp=sht["B23"]#ac platform
    return ks
end=#
