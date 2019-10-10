#=
This file contains functions that manipulate the equipment data
Sections:
-Cables
-transformers
-Per unit
=#
################################################################################
######################### Cables ###############################################
################################################################################
#loads values into the end of an array.
function eqpF_cbls_caps(cbls,km)
    for i in cbls
        push!(i,eqpF_km_cap(km,i[1],i[4],i[5]))
    end
    return cbls
end

#Calculates the new hvac cable capacity after 50-50 compensation at distance km.
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

#loads a value into the last position of each array in an aray of arrays
function eqpF_pushArray(eqp,array)
    for i in array
        push!(eqp,i)
    end
    return eqp
end

#Logic function that gets the cable data of appropriae voltage level
function eqpF_cbl_opt(kv,km)
    if kv==33.0
        opt=eqpD_33cbl_opt(km)
    elseif kv==66.0
        opt=eqpD_66cbl_opt(km)
    elseif kv==132.0
        opt=eqpD_132cbl_opt(km)
    elseif kv==220.0
        opt=eqpD_220cbl_opt(km)
    elseif kv==400.0
        opt=eqpD_400cbl_opt(km)
    elseif kv==150.0
        opt=eqpD_150cbl_opt()
    elseif kv==300.0
        opt=eqpD_300cbl_opt()
    else
        error("No matching cables for specified KV.")
    end
    return opt
end

#Fills in the physical data of a cable into the cable structure
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

#Selects sets of cables that satisfy ampacity requirements given by limits
function eqpF_cbl_sel(cbls,S,l)
    cbls_2use=Array{cbl,1}()
#Get limits and max cables possible in parallel - specified in eqp_data.jl
    lims=eqpD_eqp_lims()
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
#return cable inductive reactance
function eqpF_xl(l)
    xl=2*pi*eqpD_freq()*l
    return xl
end

#return cable capacitive reactance
function eqpF_yc(c)
    yc=2*pi*eqpD_freq()*c
    return yc
end

################################################################################
######################### Transformers #########################################
################################################################################
#built chosen sizes into transformer structured array
function eqpF_xfo_struct(s,num)
    xfm=xfo()
    xfm.mva=s
    xfm.num=num
    xfm.eta=eqpD_xEFF()
    eqpD_xfo_fail(xfm)#Set failure data
    return xfm
end

#Selects sets of transformers that satisfy power requirements given limits
function eqpF_xfo_sel(xfos,S)
    xfms_2use=Array{xfo,1}()
    lims=eqpD_eqp_lims()#Get limits and max xfos possible in parallel - specified in eqp_data.jl
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

################################################################################
################################## HVDC Converter ##############################
################################################################################
function eqpF_cnv_sel(S)
    cnv_2use=Array{xfo,1}()
    mxCnv=eqpD_cnv_mx()
    mnNum=trunc(Int,ceil(S/mxCnv))
    #Smallest possible number of converters for transmisison requirements
    cnv1=xfo()
    cnv1.mva=S/mnNum
    cnv1.num=mnNum
    cnv1.eta=eqpD_recEta()
    eqpD_cnv_fail(cnv1)
    #Add another converter
    cnv2=xfo()
    cnv2.mva=S/(mnNum+1)
    cnv2.num=mnNum+1
    cnv2.eta=eqpD_recEta()
    eqpD_cnv_fail(cnv2)
    #Add another converter
    cnv3=xfo()
    cnv3.mva=S/(mnNum+2)
    cnv3.num=mnNum+2
    cnv3.eta=eqpD_recEta()
    eqpD_cnv_fail(cnv3)
    push!(cnv_2use,cnv1)
    push!(cnv_2use,cnv2)
    push!(cnv_2use,cnv3)
    return cnv_2use
end
################################################################################
######################### Per Unit #############################################
################################################################################
#sets PU values
function eqpF_pu()
    va=lod_cnceMva()*10^6
    v=lod_pccKv()*10^3
    i=(va)/(sqrt(3)*v)
    z=v^2/va
    y=1/z
    return [va,v,i,z,y]
end

#change impedance to PU values
function eqpF_puImped(l,cb)
    cb.elec.ohm=((cb.elec.ohm*l)/cb.num)/eqpF_pu()[4]
    cb.elec.xl=((cb.elec.xl*l)/cb.num)/eqpF_pu()[4]
    cb.elec.yc=cb.elec.yc*l*cb.num*eqpF_pu()[4]
end

#Changes the base of PU value
function eqpF_puChgBs(Sn,z)
    z=z*100*10^6/Sn
    return z
end

#converts transformer and cable values to PU
function eqpF_xcXrlPu(l,xfm)
	xfm.elec.ohm=(xfm.elec.ohm/xfm.num)/eqpF_pu()[4]
	xfm.elec.xl=(xfm.elec.xl/xfm.num)/eqpF_pu()[4]
end


#combines and summarizes PU impedances
function eqpF_sumryPu(l,kv,cbx)
    eqpD_xfoXR(kv,cbx.xfm)
    eqpF_xcXrlPu(l,cbx.xfm)
    eqpF_puImped(l,cbx.cable)
    cbx.cable.elec.ohm=cbx.cable.elec.ohm+cbx.xfm.elec.ohm
    cbx.cable.elec.xl=cbx.cable.elec.xl+cbx.xfm.elec.xl
end
