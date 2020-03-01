include("wind/functions.jl")#
include("equipment/functions.jl")#
include("eens/functions.jl")#
#**
################################################################################
#################################### export functions ##########################
#given a MV a length of cable, power, wind calcs cost
function cst_mvSingleCble(cabl,S,ks,wp)
    cabl.costs.cbc=cstF_cbl_cpx(cabl)#capex
    cabl.costs.qc=0.0#assuming wind turbines can meet all q
    cabl.costs.rlc=cstF_acCbl_rlc(cabl,S,ks,wp)#cost of losses
    cabl.costs.cm=cstF_eqp_cm(cabl,ks)#corrective maintenance
    cabl.costs.eens=eensF_eqp_eens(cabl,S,ks,wp)#eens calculation
    cabl.costs.ttl=cstF_cbl_sum(cabl.costs)#totals the cable cost
    return cabl
end

#given a HV a length of cable, power, wind calcs cost
function cst_hvSingleCble(cabl,S,ks,wp)
    cabl.costs.cbc=cstF_cbl_cpx(cabl)#capex
    cabl.costs.qc=cstF_cbl_qo2o(cabl,ks)#cost of compensastion
    cabl.costs.rlc=cstF_acCbl_rlc(cabl,S,ks,wp)#cost of losses
    cabl.costs.cm=cstF_eqp_cm(cabl,ks)#corrective maintenance
    cabl.costs.eens=eensF_eqp_eens(cabl,S,ks,wp)#eens calculation
    cabl.costs.ttl=cstF_cbl_sum(cabl.costs)#totals the cable cost
    return cabl
end

#Finds max range of MV vs HV cable
function cstF_mVrng(rd,S,wp,ks,sys,ed)
    l=rd
    cst_mv=0
    cst_hv=Inf
    while cst_mv < cst_hv
        l=l+sys.mvCl
        cst_mv=cstF_MvCbl(l,S,wp,ks,ed.cbls66kV).costs.ttl
        cst_hv=cstF_HvCblo2o(l-sys.mvCl,S,wp,ks,ed.cbls220kV).costs.ttl+cstF_MvCbl(sys.mvCl,S,wp,ks,ed.cbls66kV).costs.ttl+ks.FC_bld
    end
    return l-sys.mvCl
end

#Finds which MV cable and returns it with cost
function cstF_MvCbl(l,S,wp,ks,cd66)
    cb=cbl()#create 1 object of type cbl_costs
    cb.costs.ttl=Inf#Initialize to very high total for comparison
    cbls_all=eqpF_cbls_caps(deepcopy(cd66),l)
    cbls_2use=eqpF_cbl_sel(cbls_all,S,l)#Selects 1 to ...(adjust in data) of the cables in parallel appropriate for required capacity
    for value in cbls_2use
        value=cst_mvSingleCble(value,S,ks,wp)
        if value.costs.ttl<cb.costs.ttl
            cb=deepcopy(value)#store lowest cost option
        end
    end
    return cb#return optimal cable object
end



#Finds which HV cable from oss 2 oss and returns it with cost
function cstF_HvCblo2o(l,S,wp,ks,cd220)
    cb=cbl()#create 1 object of type cbl_costs
    cb.costs.ttl=Inf#Initialize to very high total for comparison
    cbls_all=eqpF_cbls_caps(deepcopy(cd220),l)#returns all base data available for kv cables
    cbls_2use=eqpF_cbl_sel(cbls_all,S,l)#Selects 1 to ...(adjust in data) of the cables in parallel appropriate for required capacity
    for value in cbls_2use
        value=cst_hvSingleCble(value,S,ks,wp)
        if value.costs.ttl<cb.costs.ttl
            cb=deepcopy(value)#store lowest cost option
        end
    end
    return cb#return optimal cable object
end

################################################################################
############################## Cable Costs #####################################
#CAPEX of cable **
function cstF_cbl_cpx(cbl)
    cbc=cbl.length*cbl.num*cbl.cost
    return cbc
end

#AC cable loss cost **
function cstF_acCbl_rlc(cbl,s,ks,wp)
    eta=eqpD_xEFF()
    A=s*eta#eta is the efficiencu of a feeder transformer, s power to transmit
    B=cbl.elec.volt*sqrt(3)
    I=A/B#cable current
    R=(cbl.length*cbl.elec.ohm)/(cbl.num)#cable resistance
    #I^2R losses times cost factors
    #delta is related to wind profile, T_op lifetime hours and E_op cost of energy
    rlc=I^2*R*ks.T_op*ks.E_op*wp.delta
    return rlc
end

#corrective maintenance of equipment
#eq is the equipment, k are the cost factors **
function cstF_eqp_cm(eq,k)
    A=(eq.num*eq.reliability.mc)
    B=(1/eq.reliability.fr)
    C=(eq.reliability.mttr*30.417*24.0)/8760.0#changes time base of mttr of equipment: mnth-hrs/yr-hrs
    cm=k.cf*(A/(B+C))
    return cm
end

#sums all cable costs and returns the total**
function cstF_cbl_sum(res)
    ttl=res.rlc+res.qc+res.cbc+res.cm+res.eens
    return ttl
end

#AC compensation cost **
function cstF_cbl_qo2o(cbl,ks)
#div sets compensation to 50-50 split
    f=eqpD_freq()
    div=0.5
    A=cbl.elec.farrad*cbl.length*cbl.num
    Q=2*pi*f*cbl.elec.volt^2*A
    Q_oss=Q*div
    Q_pcc=Q*(1-div)
    qc=ks.Qc_oss*Q_oss+ks.Qc_oss*Q_pcc
    return qc
end

################################################################################
################### Following should be brought out to excel file ##############

#the efficiency of transformers **
function eqpD_xEFF()
    eta=0.994
    return eta
end

#set the system AC frequency **
function eqpD_freq()
    return 50.0
end

#Sets the limits that cables will be sized as a % of OWPP capacity **
function eqpD_eqp_lims(S)
        range=[0.99,1.6]
    return range
end

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
