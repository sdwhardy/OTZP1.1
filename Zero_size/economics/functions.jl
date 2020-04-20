include("wind/functions.jl")#
include("equipment/functions.jl")#
include("eens/functions.jl")#
#**
################################################################################
#################################### exported functions ########################

############################# returns cable and transformer/con... #############
#**
#Find HVAC or HVDC PCC export cable and transformer/converter (DC side not yet done) - no known voltage level
function cstF_HvCblallKvo2p(l,S,wp,ks,pcc,eqd)
    x0=xfo()
    x0.costs.ttl=0
    xPcc=cstF_xfo_pcc(S,wp,ks)
    x0.hv=pcc.kV
    xPcc.hv=pcc.kV
    c132=cstF_HvCblo2p(l,S,wp,ks,eqd.cbls132kV)
    if (132 == pcc.kV)
        x132=deepcopy(x0)
    else
        x132=deepcopy(xPcc)
    end
    x132.lv=132
    c220=cstF_HvCblo2p(l,S,wp,ks,eqd.cbls220kV)
    if (220 == pcc.kV)
        x220=deepcopy(x0)
    else
        x220=deepcopy(xPcc)
    end
    x220.lv=220
    c400=cstF_HvCblo2p(l,S,wp,ks,eqd.cbls400kV)
    if (400 == pcc.kV)
        x400=deepcopy(x0)
    else
        x400=deepcopy(xPcc)
    end
    x400.lv=400
    xfs=[x132,x220,x400]
    cbls=[c132,c220,c400]
    low_cost=findmin([c132.costs.ttl+x132.costs.ttl,c220.costs.ttl+x220.costs.ttl,c400.costs.ttl+x400.costs.ttl])[2]
    return cbls[low_cost],xfs[low_cost]
end

#** - Finds owpp to OSS cables/transformers - does check for fixed construction cost of OSS
function cstF_MvHvCbloss(l,S,wp,ks,kv,sys,mv_mx,eqd)
    if (l<=mv_mx)
        #MV
        cmv=cstF_MvCbl(l,S,wp,ks,eqd.cbls66kV)
        xfm0=cstF_xfo_oss(S,wp,ks)
        xfm0.hv=kv
        xfm0.lv=66
        return [cmv],[xfm0]
    else
        #HV
        cmv_05=cstF_MvCbl(sys.mvCl,S,wp,ks,eqd.cbls66kV)
        chv,xhv=cstF_HvCblallKvo2o(l-sys.mvCl,S,wp,ks,kv,eqd)
        xhv[1].lv=cmv_05.elec.volt
        xhv[2].lv=chv.elec.volt
        return [cmv_05,chv],[xhv[1],xhv[2]]
    end
end

function cstF_MvHvCblpcc(l,S,wp,ks,pcc,sys,mv_mx,eqd)
    if (l<=mv_mx)
        #MV
        cmv=cstF_MvCbl(l,S,wp,ks,eqd.cbls66kV)
        xfm0=cstF_xfo_pcc(S,wp,ks)
        xfm0.hv=pcc.kV
        xfm0.lv=66
        return [cmv,xfm0]
    else
        #HV
        cmv_05=cstF_MvCbl(sys.mvCl,S,wp,ks,eqd.cbls66kV)
        chv,xhv=cstF_HvCblallKvo2p(l-sys.mvCl,S,wp,ks,pcc,eqd)
        xhv.lv=cmv_05.elec.volt
        xfm0=cstF_xfo_oss(S,wp,ks)
        return [cmv_05,chv,xfm0,xhv]
    end
end
#**
#returns OSS to OSS HV cable and required trandformers when voltage is unknown
function cstF_HvCblallKvo2o(l,S,wp,ks,oss_kv,eqd)
    #set transformers
    xfm0=xfo()
    xfm0.costs.ttl=0
    xfm0.mva=0
    xfm0.num=0
    xfm0.hv=oss_kv
    xfm1=cstF_xfo_oss(S,wp,ks)
    xfm1.hv=oss_kv

    if (oss_kv==132)
        x132=deepcopy(xfm0)
    else
        x132=deepcopy(xfm1)
    end
    x132.lv=132
    if (oss_kv==220)
        x220=deepcopy(xfm0)
    else
        x220=deepcopy(xfm1)
    end
    x220.lv=220
    if (oss_kv==400)
        x400=deepcopy(xfm0)
    else
        x400=deepcopy(xfm1)
    end
    x400.lv=400
    xs=[x132,x220,x400]

    #set cable
    c132=cstF_HvCblo2o(l,S,wp,ks,eqd.cbls132kV)
    c220=cstF_HvCblo2o(l,S,wp,ks,eqd.cbls220kV)
    c400=cstF_HvCblo2o(l,S,wp,ks,eqd.cbls400kV)
    cbls=[c132,c220,c400]
    low_cost=findmin([c132.costs.ttl+x132.costs.ttl,c220.costs.ttl+x220.costs.ttl,c400.costs.ttl+x400.costs.ttl])[2]
    return cbls[low_cost],[xfm1,xs[low_cost]]
end

################################# Transformers #################################
#**
#**
function opt_mogXfmrs(mog,pMv,pHv,p132,p220,p400,wMv,wHv,w132,w220,w400,ks,kV)
    wo=wind()
    wo.pu=zeros(Float32,8759)
    wo.ce=zeros(Float32,8759)
    wo.delta=0
    wo.lf=0
    if length(pMv) !=0
        push!(mog.xfmrs,cstF_xfo_oss(sum(pMv),opt_Wsum(deepcopy(wo),wMv),ks,66,kV))
    end
    if length(p132) !=0
        push!(mog.xfmrs,cstF_xfo_oss(sum(p132),opt_Wsum(deepcopy(wo),w132),ks,132,kV))
    end
    if length(p220) !=0
        push!(mog.xfmrs,cstF_xfo_oss(sum(p220),opt_Wsum(deepcopy(wo),w220),ks,220,kV))
    end
    if length(p400) !=0
        push!(mog.xfmrs,cstF_xfo_oss(sum(p400),opt_Wsum(deepcopy(wo),w400),ks,400,kV))
    end
    return mog
end

#finds best transformer from all options at OSS
function cstF_xfo_oss(S,wp::wind,ks::cstS_ks,lv=66.0,hv=66.0)
    xfm=xfo()#create transformer object
    xfm.costs.ttl=Inf
    xfos_all=eqpD_xfo_opt()#get all available xformer sizes
    xfos_2use=eqpF_xfo_sel(xfos_all,S)#select combinations to calculate
    for valu in xfos_2use
        valu.lv=lv
        valu.hv=hv
        valu=cst_ossXfmr(valu,S,ks,wp)
        #store lowest cost option
        if valu.costs.ttl<xfm.costs.ttl
            xfm=deepcopy(valu)
        end
    end
    return xfm
end

#finds best transformer from all options transformer at PCC
function cstF_xfo_pcc(S,wp,ks)
    xfm=xfo()#create transformer object
    xfm.costs.ttl=Inf
    xfos_all=eqpD_xfo_opt()#get all available xformer sizes
    xfos_2use=eqpF_xfo_sel(xfos_all,S)#select combinations to calculate

    for value in xfos_2use
        value=cst_pccXfmr(value,S,ks,wp)
        #store lowest cost option
        if (value.costs.ttl<xfm.costs.ttl)
            xfm=deepcopy(value)
        end
    end
    return xfm
end

#OSS transformer set cost
function cst_ossXfmr(xfmr,S,ks,wp)
    xfmr.costs.cpx=cstF_oss_cpx(xfmr,ks)#capex oss
    xfmr.costs.tlc=cstF_xfo_tlc(xfmr,S,ks,wp)#cost of losses
    xfmr.costs.cm=cstF_eqp_cm(xfmr,ks)#corrective maintenance
    xfmr.costs.eens=eensF_eqp_eens(xfmr,S,ks,wp)#eens calculation
    xfmr.costs.ttl=cstF_xfo_sum(xfmr.costs)#totals the xfo cost
    return xfmr
end

#PCC transformer set cost
function cst_pccXfmr(xfmr,S,ks,wp)
    xfmr.costs.cpx=cstF_pcc_cpx(xfmr)#capex oss
    xfmr.costs.tlc=cstF_xfo_tlc(xfmr,S,ks,wp)#cost of losses
    xfmr.costs.cm=cstF_eqp_cm(xfmr,ks)#corrective maintenance
    xfmr.costs.eens=eensF_eqp_eens(xfmr,S,ks,wp)#eens calculation
    xfmr.costs.ttl=cstF_xfo_sum(xfmr.costs)#totals the xfo cost
    return xfmr
end
#################################### mv Cables Only ############################
#Finds max range of MV vs HV cable
function cstF_mVrng(rd,S,wp,ks,sys,ed)
    l=rd
    cst_mv=0
    cst_hv=Inf
    if (sys.mvCl>=0.01)
        increment=sys.mvCl
    else
        increment=0.01
    end
    while cst_mv < cst_hv
        l=l+increment
        cst_mv=cstF_MvCbl(l,S,wp,ks,ed.cbls66kV).costs.ttl
        cst_hv=cstF_HvCblo2o(l-increment,S,wp,ks,ed.cbls220kV).costs.ttl+cstF_MvCbl(increment,S,wp,ks,ed.cbls66kV).costs.ttl+ks.FC_bld
    end
    return l-increment
end

#MV cables
#Finds which OWPP MV cable from ALL OPTIONS
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

#given a MV a length of cable, power, wind calcs cost
function cst_mvSingleCble(cabl,S,ks,wp)
    cabl.costs.cpx=cstF_cbl_cpx(cabl)#capex
    cabl.costs.qc=0.0#assuming wind turbines can meet all q
    cabl.costs.rlc=cstF_acCbl_rlc(cabl,S,ks,wp)#cost of losses
    cabl.costs.cm=cstF_eqp_cm(cabl,ks)#corrective maintenance
    cabl.costs.eens=eensF_eqp_eens(cabl,S,ks,wp)#eens calculation
    cabl.costs.ttl=cstF_cbl_sum(cabl.costs)#totals the cable cost
    return cabl
end

#################################### hv Cables Only ############################
#Finds which HV cable from oss 2 oss and returns it with cost  oss 2 oss
function cstF_HvCblo2o(l,S,wp,ks,cd220)
    cb=cbl()#create 1 object of type cbl_costs
    cb.costs.ttl=Inf#Initialize to very high total for comparison
    cbls_all=eqpF_cbls_caps(deepcopy(cd220),l)#returns all base data available for kv cables
    cbls_2use=eqpF_cbl_sel(cbls_all,S,l)#Selects 1 to ...(adjust in data) of the cables in parallel appropriate for required capacity
    for value in cbls_2use
        value=cst_hvSingleCbleo2o(value,S,ks,wp)
        if value.costs.ttl<cb.costs.ttl
            cb=deepcopy(value)#store lowest cost option
        end
    end
    return cb#return optimal cable object
end

#Finds which HV cable from oss 2 pcc and returns it with cost  oss 2 pcc
function cstF_HvCblo2p(l,S,wp,ks,cd220)
    cb=cbl()#create 1 object of type cbl_costs
    cb.costs.ttl=Inf#Initialize to very high total for comparison
    cbls_all=eqpF_cbls_caps(deepcopy(cd220),l)#returns all base data available for kv cables
    cbls_2use=eqpF_cbl_sel(cbls_all,S,l)#Selects 1 to ...(adjust in data) of the cables in parallel appropriate for required capacity
    for value in cbls_2use
        value=cst_hvSingleCbleo2p(value,S,ks,wp)
        if value.costs.ttl<cb.costs.ttl
            cb=deepcopy(value)#store lowest cost option
        end
    end
    return cb#return optimal cable object
end

#HV cables
#given a HV a length of cable, power, wind calcs cost  oss 2 oss
function cst_hvSingleCbleo2o(cabl,S,ks,wp)
    cabl.costs.cpx=cstF_cbl_cpx(cabl)#capex
    cabl.costs.qc=cstF_cbl_qo2o(cabl,ks)#cost of compensastion
    cabl.costs.rlc=cstF_acCbl_rlc(cabl,S,ks,wp)#cost of losses
    cabl.costs.cm=cstF_eqp_cm(cabl,ks)#corrective maintenance
    cabl.costs.eens=eensF_eqp_eens(cabl,S,ks,wp)#eens calculation
    cabl.costs.ttl=cstF_cbl_sum(cabl.costs)#totals the cable cost
    return cabl
end

#given a HV a length of cable, power, wind calcs cost oss 2 pcc
function cst_hvSingleCbleo2p(cabl,S,ks,wp)
    cabl.costs.cpx=cstF_cbl_cpx(cabl)#capex
    cabl.costs.qc=cstF_cbl_qo2p(cabl,ks)#cost of compensastion
    cabl.costs.rlc=cstF_acCbl_rlc(cabl,S,ks,wp)#cost of losses
    cabl.costs.cm=cstF_eqp_cm(cabl,ks)#corrective maintenance
    cabl.costs.eens=eensF_eqp_eens(cabl,S,ks,wp)#eens calculation
    cabl.costs.ttl=cstF_cbl_sum(cabl.costs)#totals the cable cost
    return cabl
end

################################################################################
############################## Cable Costs #####################################
#CAPEX of cable **
function cstF_cbl_cpx(cbl)
    cpx=cbl.length*cbl.num*cbl.cost
    return cpx
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
    ttl=res.rlc+res.qc+res.cpx+res.cm+res.eens
    return ttl
end

#AC compensation cost oss 2 oss**
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

#AC compensation cost oss 2 pcc**
function cstF_cbl_qo2p(cbl,ks)
#div sets compensation to 50-50 split
    f=eqpD_freq()
    div=0.5
    A=cbl.elec.farrad*cbl.length*cbl.num
    Q=2*pi*f*cbl.elec.volt^2*A
    Q_oss=Q*div
    Q_pcc=Q*(1-div)
    qc=ks.Qc_oss*Q_oss+ks.Qc_pcc*Q_pcc
    return qc
end

################################################################################
############################### transformer cost ###############################
#OSS platform CAPEX
#excluding fixed platform cost**
function cstF_oss_cpx(x,k)
    A=(1+k.dc*(x.num-2))
    B=(k.f_ct+k.p_ct)
    #B=0.075
    cpx=A*B*x.num*x.mva
    return cpx
end
#PCC transformer **
function cstF_pcc_cpx(x)
    cpx=0.02621*(x.mva*x.num)^0.7513
    return cpx
end
#xfm Losses Calculation**
function cstF_xfo_tlc(xfo,S,ks,wp)
    pf=eqpD_pf()
    tlc=S*pf*(1-xfo.eta)*ks.T_op*ks.E_op*wp.delta
    return tlc
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

#sums all transformer/converter costs and returns the total**
function cstF_xfo_sum(res)
    ttl=res.cpx+res.tlc+res.cm+res.eens
    return ttl
end

################################################################################
#returns set of kV cables
function cstF_kvCblSet(kV,eqd)
    if (kV==33)
        cdkV=eqd.cbls33kV
    elseif (kV==66)
        cdkV=eqd.cbls66kV
    elseif (kV==132)
        cdkV=eqd.cbls132kV
    elseif (kV==220)
        cdkV=eqd.cbls220kV
    elseif (kV==400)
        cdkV=eqd.cbls400kV
    else
        println("Error: no kV cable level match!!!!!!!!!!!!!!")
    end
    return cdkV
end

#checks 1 size down in MV cable
function cstF_MvCbl_nextSizeDown(l,S,kv,wp,ks,cabl,nm,eqd)
    cb=cbl()#create 1 object of type cbl_costs
    cb.costs.ttl=Inf#Initialize to very high total for comparison
    cdkV=cstF_kvCblSet(kv,eqd)
    cbls_2use=eqpF_nextSizeDown(S,l,cabl,cdkV,nm)
    for value in cbls_2use
        value=cst_mvSingleCble(value,S,ks,wp)
        if value.costs.ttl<cb.costs.ttl
            cb=deepcopy(value)#store lowest cost option
        end
    end
    return cb#return optimal cable object
end

#checks 1 size up in MV cable
function cstF_MvCbl_nextSizeUp(l,S,kv,wp,ks,cabl,nm,eqd)
    cb=cbl()#create 1 object of type cbl_costs
    cb.costs.ttl=Inf#Initialize to very high total for comparison
    cdkV=cstF_kvCblSet(kv,eqd)
    cbls_2use=eqpF_nextSizeUp(S,l,cabl,cdkV,nm)#Selects 1 to ...(adjust in data) of the cables in parallel appropriate for required capacity
    for value in cbls_2use
        value=cst_mvSingleCble(value,S,ks,wp)
        if value.costs.ttl<cb.costs.ttl
            cb=deepcopy(value)#store lowest cost option
        end
    end
    return cb#return optimal cable object
end

#checks 1 size down in HV cable
function cstF_nextSizeDown(l,S,kv,wp,ks,cabl,nm,eqd)
    cb=cbl()#create 1 object of type cbl_costs
    cb.costs.ttl=Inf#Initialize to very high total for comparison
    cdkV=cstF_kvCblSet(kv,eqd)
    cbls_2use=eqpF_nextSizeDown(S,l,cabl,cdkV,nm)#Selects 1 to ...(adjust in data) of the cables in parallel appropriate for required capacity
    for value in cbls_2use
        value=cst_hvSingleCbleo2o(value,S,ks,wp)
        if value.costs.ttl<cb.costs.ttl
            cb=deepcopy(value)#store lowest cost option
        end
    end
    return cb#return optimal cable object
end

#checks 1 size up in HV cable
function cstF_nextSizeUp(l,S,kv,wp,ks,cabl,nm,eqd)
    cb=cbl()#create 1 object of type cbl_costs
    cb.costs.ttl=Inf#Initialize to very high total for comparison
    cdkV=cstF_kvCblSet(kv,eqd)
    cbls_2use=eqpF_nextSizeUp(S,l,cabl,cdkV,nm)#Selects 1 to ...(adjust in data) of the cables in parallel appropriate for required capacity
    for value in cbls_2use
        value=cst_hvSingleCbleo2o(value,S,ks,wp)
        if value.costs.ttl<cb.costs.ttl
            cb=deepcopy(value)#store lowest cost option
        end
    end
    return cb#return optimal cable object
end

#checks 1 size down in PCC cable
function cstF_nextSizeDownPcc(l,S,kv,wp,ks,cabl,nm,eqd)
    cb=cbl()#create 1 object of type cbl_costs
    cb.costs.ttl=Inf#Initialize to very high total for comparison
    cdkV=cstF_kvCblSet(kv,eqd)
    cbls_2use=eqpF_nextSizeDown(S,l,cabl,cdkV,nm)#Selects 1 to ...(adjust in data) of the cables in parallel appropriate for required capacity
    for value in cbls_2use
        value=cst_hvSingleCbleo2p(value,S,ks,wp)
        if value.costs.ttl<cb.costs.ttl
            cb=deepcopy(value)#store lowest cost option
        end
    end
    return cb#return optimal cable object
end

#checks 1 size up in PCC cable
function cstF_nextSizeUpPcc(l,S,kv,wp,ks,cabl,nm,eqd)
    cb=cbl()#create 1 object of type cbl_costs
    cb.costs.ttl=Inf#Initialize to very high total for comparison
    cdkV=cstF_kvCblSet(kv,eqd)
    cbls_2use=eqpF_nextSizeUp(S,l,cabl,cdkV,nm)#Selects 1 to ...(adjust in data) of the cables in parallel appropriate for required capacity
    for value in cbls_2use
        value=cst_hvSingleCbleo2p(value,S,ks,wp)
        if value.costs.ttl<cb.costs.ttl
            cb=deepcopy(value)#store lowest cost option
        end
    end
    return cb#return optimal cable object
end
