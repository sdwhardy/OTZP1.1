################################################################################
############################## Cable Structure #################################
#loads values into the end of an array. **
function eqpF_cbls_caps(cbls,km)
    for i in cbls
        push!(i,eqpF_km_cap(km,i[1],i[4],i[5]))
    end
    return cbls
end

#Calculates the new hvac cable capacity after 50-50 compensation at distance km. **
function eqpF_km_cap(l,v,q,a)
#get system frequency
    f=eqpD_freq()
#Calculates the square of new hvac cable capacity after 50-50 compensation at distance km.
    mva=(sqrt(3)*v*10^3*a/10^6)^2-((0.5*((v*10^3)^2*2*pi*f*l*q*10^-9))/10^6)^2
#takes square root if positive returns zero if negative
    if mva>=0
        mva=sqrt(mva)
    else
        mva=0.0
    end
 return mva
end

#Fills in the physical data of a cable into the cable structure **
function eqpF_cbl_struct(cb,km,num)
    cbl_data=cbl()
    cbl_data.elec.volt=cb[1]
    cbl_data.size=cb[2]
    cbl_data.elec.ohm=cb[3]*10^-3
    cbl_data.elec.farrad=cb[4]*10^-9
    cbl_data.elec.amp=cb[5]
    cbl_data.cost=cb[6]*10^-3
    cbl_data.length=km
    cbl_data.elec.henry=cb[7]*10^-3
    cbl_data.elec.xl=eqpF_xl(cbl_data.elec.henry)
    cbl_data.elec.yc=eqpF_yc(cbl_data.elec.farrad)
    cbl_data.mva=cb[8]
    cbl_data.num=num
    eqpD_cbl_fail(cbl_data)#Set failure data
    return cbl_data
end

#Selects sets of cables that satisfy ampacity requirements given by limits **
function eqpF_cbl_sel(cbls,S,l)
    cbls_2use=Array{cbl,1}()
#Get limits and max cables possible in parallel - specified in eqp_data.jl
    lims=eqpD_eqp_lims(S)
    parCmax=eqpD_MAXcbls(cbls[1][1])
    for i in cbls
        for j=1:parCmax
            if ((j*i[8])>lims[1]*S && (j*i[8])<lims[2]*S)
                push!(cbls_2use,eqpF_cbl_struct(i,l,j))
            end
        end
    end
    return cbls_2use
end

#return cable inductive reactance **
function eqpF_xl(l)
    xl=2*pi*eqpD_freq()*l
    return xl
end

#return cable capacitive reactance **
function eqpF_yc(c)
    yc=2*pi*eqpD_freq()*c
    return yc
end


################################################################################
################### Following should be brought out to excel file ##############
#failure data for cables **
function eqpD_cbl_fail(cbl)
    cbl.reliability.fr=(0.08/100)*cbl.length#/yr/100km
    cbl.reliability.mttr=2.0#/yr/100km
    cbl.reliability.mc=0.56
    return nothing
end

#Set maximum of cables possible in parallel **
function eqpD_MAXcbls(kv)
    if kv == 33 || kv == 66
        pll=12
    else
        pll=15
    end
    return pll
end

################################################################################
######################### Transformers #########################################
################################################################################
#Selects sets of transformers that satisfy power requirements given limits **
function eqpF_xfo_sel(xfos,S)
    xfms_2use=Array{xfo,1}()
    lims=eqpD_eqp_lims(S)#Get limits and max xfos possible in parallel - specified in eqp_data.jl
    parXmax=eqpD_MAXxfos()
    for i in xfos
        for j=1:parXmax
            if ((j*i)>lims[1]*S && (j*i)<lims[2]*S)
                push!(xfms_2use,eqpF_xfo_struct(i,j))
            end
        end
    end
    return xfms_2use
end

#built chosen sizes into transformer structured array **
function eqpF_xfo_struct(s,num)
    xfm=xfo()
    xfm.mva=s
    xfm.num=num
    xfm.eta=eqpD_xEFF()
    eqpD_xfo_fail(xfm)#Set failure data
    return xfm
end
################################################################################
################### Following should be brought out to excel file ##############
#failure data for transformers **
function eqpD_xfo_fail(x)
    x.reliability.fr=0.03#/yr
    x.reliability.mttr=6.0#month
    x.reliability.mc=2.8#
    return nothing
end

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

################################################################################
################################ General #######################################
################################################################################
################### Following should be brought out to excel file ##############

#set the system AC frequency **
function eqpD_freq()
    return 50.0
end

#Sets the limits that cables will be sized as a % of OWPP capacity **
function eqpD_eqp_lims(S)
        range=[0.99,1.6]
    return range
end

#sets owpp power factor **
function eqpD_pf()
    return 1.0
end

################################################################################

function eqpF_nextSizeDown(S,l,sze,cdkV,nm)
    cbls_all=eqpF_cbls_caps(deepcopy(cdkV),l)
    cbls_2use=eqpF_cbl_sel(cbls_all,S,l)#Selects 1 to ...(adjust in data) of the cables in parallel appropriate for required capacity
    cbls_2return=Array{cbl,1}()
    cbl0=cbl()
    cbl1=cbl()
    for cb in cbls_2use
        if cb.size==sze && cb.num==nm
            cbl0=deepcopy(cb)
            break
        end
    end

    cbl1.mva=0
    cbl1.num=0

    for cb in cbls_2use
        if (cb.mva*cb.num<cbl0.mva*cbl0.num && cb.mva*cb.num>cbl1.mva*cbl1.num)
            cbl1=deepcopy(cb)
        end
    end
    push!(cbls_2return,cbl0)
    push!(cbls_2return,cbl1)
    return cbls_2return
end

#sze=cabl
function eqpF_nextSizeUp(S,l,sze,cdkV,nm)
    cbls_all=eqpF_cbls_caps(deepcopy(cdkV),l)
    cbls_2use=eqpF_cbl_sel(cbls_all,S,l)#Selects 1 to ...(adjust in data) of the cables in parallel appropriate for required capacity
    cbls_2return=Array{cbl,1}()
    cbl0=cbl()
    cbl1=cbl()
    #cb=cbls_2use[length(cbls_2use)-3]
    for cb in cbls_2use
        if (cb.size==sze && cb.num==nm)
            cbl0=deepcopy(cb)
            break
        end
    end

    cbl1.mva=Inf
    cbl1.num=Inf

    for cb in cbls_2use
        if (cb.mva*cb.num>cbl0.mva*cbl0.num && cb.mva*cb.num<cbl1.mva*cbl1.num)
            cbl1=deepcopy(cb)
        end
    end
    push!(cbls_2return,deepcopy(cbl0))
    push!(cbls_2return,deepcopy(cbl1))
    return cbls_2return
end
