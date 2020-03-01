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
#cbls=eqpF_cbl_opt(kv,l)
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
