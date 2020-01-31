#=
This file contains the input data for equipment used
Sections:
-cables
-transformers
-converters
-grid
=#
################################################################################
######################### Cables ###############################################
################################################################################

#Set maximum of cables possible in parallel **
function eqpD_MAXcbls(kv)
    if kv == 33 || kv == 66
        pll=12
    else
        pll=15
    end
    return pll
end

#km=50
#33kV cables **
function eqpD_33cbl_opt(km)
    p2e=cstD_xchg()#exchange rate
    #%kV,cm^2,mohms/km,nF/km,Amps,10^3 euros/km, mH, capacity at km
    a=[33, 95, 218, 0.18, 300, 434*p2e, 0.44]
    b=[33, 120, 172, 0.19, 340, 443*p2e, 0.42]
    c=[33, 150, 136, 0.21, 375, 453*p2e, 0.41]
    d=[33, 185, 110, 0.22, 420, 466*p2e, 0.39]
    e=[33, 240, 84.8, 0.24, 480, 468*p2e, 0.38]
    f=[33, 300, 67.6, 0.26, 530, 505*p2e, 0.36]
    g=[33, 400, 53.2, 0.29, 590, 531*p2e, 0.35]
    h=[33, 500, 42.8, 0.32, 655, 564*p2e, 0.34]
    i=[33, 630, 34.6, 0.35, 715, 598*p2e, 0.32]
    j=[33, 800, 28.7, 0.38, 775, 638*p2e, 0.31]
    alphbt=[a,b,c,d,e,f,g,h,i,j]
    alphbt=eqpF_cbls_caps(alphbt,km)
    return alphbt
end

#66kV cables **
function eqpD_66cbl_opt(km)
    p2e=cstD_xchg()#exchange rate
    #%kV,cm^2,mohms/km,nF/km,Amps,10^3 euros/km, mH, capacity at km
    a=[66, 95, 218, 0.17, 300, 462*p2e, 0.44]
    b=[66, 120, 172, 0.18, 340, 472*p2e, 0.43]
    c=[66, 150, 136, 0.19, 375, 482*p2e, 0.41]
    d=[66, 185, 110, 0.2, 420, 496*p2e, 0.4]
    e=[66, 240, 84.8, 0.22, 480, 517*p2e, 0.38]
    f=[66, 300, 67.6, 0.24, 530, 537*p2e, 0.37]
    g=[66, 400, 53.2, 0.26, 590, 564*p2e, 0.35]
    h=[66, 500, 42.8, 0.29, 655, 598*p2e, 0.34]
    i=[66, 630, 34.6, 0.32, 715, 634*p2e, 0.33]
    j=[66, 800, 28.7, 0.35, 775, 676*p2e, 0.32]
    k=[66, 1000, 24.5, 0.38, 825, 716*p2e, 0.31]
    alphbt=[a,b,c,d,e,f,g,h,i,j,k]
    alphbt=eqpF_cbls_caps(alphbt,km)
    return alphbt
end

#132kV cables **
function eqpD_132cbl_opt(km)
    p2e=cstD_xchg()#exchange rate
    #%kV,cm^2,mohms/km,nF/km,Amps,10^3 euros/km, capacity at km
    a=[132,185,100,165,501,424*p2e,0.47]
    b=[132,300,76.1,175,600,504*p2e,0.42]
    c=[132,400,60.6,185,677,568*p2e,0.4]
    d=[132,500,49.3,192,739,635*p2e,0.387]
    e=[132,630,39.5,209,818,685*p2e,0.372]
    f=[132,800,32.4,217,888,795*p2e,0.364]
    g=[132,1000,27.5,238,949,860*p2e,0.351]
    alphbt=[a,b,c,d,e,f,g]
    alphbt=eqpF_cbls_caps(alphbt,km)
    return alphbt
end

#220kV cables **
function eqpD_220cbl_opt(km)
    p2e=cstD_xchg()#exchange rate
    #%kV,cm^2,mohms/km,nF/km,Amps,10^3 euros/km, capacity at km
    a=[220,400,60.1,122,665,728*p2e,0.457]
    b=[220,500,48.9,136,732,815*p2e,0.437]
    c=[220,630,39.1,151,808,850*p2e,0.415]
    d=[220,800,31.9,163,879,975*p2e,0.4]
    e=[220,1000,27,177,942,1000*p2e,0.386]
    alphbt=[a,b,c,d,e]
    alphbt=eqpF_cbls_caps(alphbt,km)
    return alphbt
end

#400kV cables **
function eqpD_400cbl_opt(km)
    p2e=cstD_xchg()#exchange rate
    cores=1
    #%kV,cm^2,mohms/km,nF/km,Amps,10^3 euros/km, capacity at km
    a=[400,500,40,117,776,1239*p2e,0.589*cores]
    b=[400,630,36,125,824,1323*p2e,0.561*cores]
    c=[400,800,31.4,130,870,1400*p2e,0.54*cores]
    d=[400,1000,26.5,140,932,1550*p2e,0.52*cores]
    e=[400,1200,22.1,170,986,1700*p2e,0.49*cores]
    f=[400,1400,18.9,180,1015,1850*p2e,0.47*cores]
    g=[400,1600,16.6,190,1036,2000*p2e,0.46*cores]
    h=[400,2000,13.2,200,1078,2150*p2e,0.44*cores]
    alphbt=[a,b,c,d,e,f,g,h]
    alphbt=eqpF_cbls_caps(alphbt,km)
    return alphbt
end

#Sets the limits that cables will be sized as a % of OWPP capacity **
function eqpD_eqp_lims(S)
    #if S<=500
        range=[0.99,1.6]
    #elseif (S>350 && S<500)
        #range=[0.95,3]
    #else
        #range=[0.95,3]
    #end
    return range
end

#failure data for cables **
function eqpD_cbl_fail(cbl)
    cbl.reliability.fr=0.04#/yr/100km
    cbl.reliability.mttr=2.0#/yr/100km
    cbl.reliability.mc=0.56
    return nothing
end

################################################################################
######################### Transformers #########################################
################################################################################
#Sets all options for transformer sizes in 10MVA steps **
function  eqpD_xfo_opt()
    xfos=Array{Float32,1}()
    push!(xfos,10)
    for i=50:10:500
        push!(xfos,i)
    end
    return xfos
end

#Set maximum of transformers possible in parallel **
function eqpD_MAXxfos()
    return 5
end

#the efficiency of transformers **
function eqpD_xEFF()
    eta=0.994
    return eta
end

#failure data for transformers **
function eqpD_xfo_fail(x)
    x.reliability.fr=0.03#/yr
    x.reliability.mttr=6.0#month
    x.reliability.mc=2.8#
    return nothing
end
################################################################################
######################### Grid #################################################
################################################################################
#sets owpp power factor **
function eqpD_pf()
    return 1.0
end

#set the system AC frequency **
function eqpD_freq()
    return 50.0
end
################################################################################
######################### Converters ###########################################
################################################################################
#=function  eqpD_cnv_mx()
    mxCnv=6400.0
    return mxCnv
end

#failure data Converters
function eqpD_cnv_fail(conv)
    conv.reliability.fr=0.12#/yr
    conv.reliability.mttr=1#month
    conv.reliability.mc=0.56#
    return nothing
end

#efficiency of rectifier
function eqpD_recEta()
    eta=0.9828
    return eta
end

#Efficiency of inverter
function eqpD_invEta()
    eta=0.9819
    return eta
end=#

#=
#150kV hvdc cables
function eqpD_150cbl_opt(cbls)
    p2e=cstD_xchg()#exchange rate
    #%kV,cm^2,mohms/km,nF/km,Amps,10^3 euros/km, mH, capacity at km
    a=[150,1000, 22.4, 0.6969, 1644, 670*p2e, 0.6969,150*2*1644*10^-3]
    b=[150,1200, 19.2, 0.6969, 1791, 730*p2e, 0.6969,150*2*1791*10^-3]
    c=[150,1400, 16.5, 0.6969, 1962, 785*p2e, 0.6969,150*2*1962*10^-3]
    d=[150,1600, 14.4, 0.6969, 2123, 840*p2e, 0.6969,150*2*2123*10^-3]
    e=[150,2000, 11.5, 0.6969, 2407, 900*p2e, 0.6969,150*2*2407*10^-3]
    alphbt=[a,b,c,d,e]
    return alphbt
end

#300kV hvdc cables
function eqpD_300cbl_opt(cbls)
    p2e=cstD_xchg()#exchange rate
    #%kV,cm^2,mohms/km,nF/km,Amps,10^3 euros/km, mH, capacity at km
    a=[300,1000, 22.4, 0.6969, 1644, 855*p2e, 0.6969,300*2*1644*10^-3]
    b=[300,1200, 19.2, 0.6969, 1791, 940*p2e, 0.6969,300*2*1791*10^-3]
    c=[300,1400, 16.5, 0.6969, 1962, 1015*p2e, 0.6969,300*2*1962*10^-3]
    d=[300,1600, 14.4, 0.6969, 2123, 1090*p2e, 0.6969,300*2*2123*10^-3]
    e=[300,2000, 11.5, 0.6969, 2407, 1175*p2e, 0.6969,300*2*2407*10^-3]
    alphbt=[a,b,c,d,e]
    return alphbt
end
=#
#admitance values for infinite/DC lines
#=function eqpD_dcAdm()
        r=0.0093
        x=0.0222
        b=0.2217
    return [r,x,b]
end=#

######################################### depricated ###########################
#=
#add percent impendance of transformer
function eqpD_xfoXR(kv,x)
    #400/132 - X=8% R=0.14% on 100MVA base
    #275/132 - X=9% R=0.16% on 100MVA base #Source National grid
    if kv == 400
        X=0.08
        R=0.0014
    elseif kv == 220 || kv == 66
        X=0.09
        R=0.0016
    #assumed values based on data given. needs to be updated
    elseif kv == 132 || kv == 33
        X=0.1
        R=0.0018
    else
        error("kV doesn't match transformer for % impedance!")
    end
    mva=eqpF_pu()[1]
    x.elec.xl=eqpF_puChgBs(mva,X)
    x.elec.ohm=eqpF_puChgBs(mva,R)
end
=#
