###################################################################
mutable struct xy
      x::Float64
      y::Float64
end
xy()=xy(69.69,69.69)
###################################################################
mutable struct gps
      lat::Float32
      lng::Float32
end
gps()=gps(69.69,69.69)
###################################################################
#wind object
mutable struct wind
      pu::Array{Float32}
      ce::Array{Float32}
      delta::Float32
      lf::Float32
end
wind()=wind(Float32[],Float32[],69.69,69.69)
###################################################################
mutable struct node
      gps::gps
      xy::xy
      num::Int32
end
node()=node(gps(),xy(),69)
###########################################################################
#data structures used for equipment are specified in this file
mutable struct elec
   amp::Float32
   volt::Float32
   ohm::Float32
   farrad::Float32
   henry::Float32
   yc::Float32
   xl::Float32
end
elec()=elec(69.69,69.69,69.69,69.69,69.69,69.69,69.69)
###########################################################################
#reliability structure
mutable struct relia
   fr::Float32
   mttr::Float32
   mc::Float32
end
relia()=relia(69.69,69.69,69.69)
###################################################################
#the structure of costs for a cable
mutable struct cbl_costs
   qc::Float32
   cpx::Float32
   rlc::Float32
   cm::Float32
   eens::Float32
   ttl::Float32
end
cbl_costs()=cbl_costs(0.0,0.0,0.0,0.0,0.0,0.0)
########################################################################
#the structure of costs for a transformers
mutable struct xfo_costs
   cpx::Float32
   tlc::Float32
   cm::Float32
   eens::Float32
   ttl::Float32
end
xfo_costs()=xfo_costs(0.0,0.0,0.0,0.0,0.0)
########################################################################
#cost components and totals calculated for a OWPP object
#mutable struct results
#     cpx::Float32
#     loss::Float32
#     opex::Float32
#     ttl::Float32
#end
#results()=results(0.0,0.0,0.0,0.0)
########################################################################
#the structure for a transformer
mutable struct xfo
   mva::Float32
   num::Float32
   hv::Float32
   lv::Float32
   eta::Float32
   reliability::relia
   elec::elec
   costs::xfo_costs
end
xfo()=xfo(0.0,0.0,0.0,0.0,0.0,relia(),elec(),xfo_costs())
###################################################################
#the structure used for a cable
mutable struct cbl
   mva::Float32
   length::Float32
   pth::Array{node}
   size::Float32
   num::Float32
   cost::Float32
   reliability::relia
   elec::elec
   costs::cbl_costs
end
cbl()=cbl(0.0,0.0,node[],0.0,0.0,0.0,relia(),elec(),cbl_costs())
###################################################################
mutable struct bus
      mva::Float32
      wnd::wind
      name::String
      node::node
      zone::Float32
      mv_zone::Float32
      kV::Int32
      num::Int32
      xfmrs::Array{xfo}
      base_cost::Float32
end
bus()=bus(69.69,wind(),"sixty-nine",node(),69.69,69.69,69,69,xfo[],10)
####################################################################
#an object that contains all cost factors used in the calculations
mutable struct cstS_ks
   FC_ac::Float32
   FC_dc::Float32
   dc::Float32
   f_ct::Float32
   p_ct::Float32
   c_ct::Float32
   Qc_oss::Float32
   Qc_pcc::Float32
   life::Float32
   T_op::Float32
   E_op::Float32
   cf::Float32
   FC_bld::Float32
   p2e::Float32
end
cstS_ks()=cstS_ks(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
########################################################################
mutable struct system
      nogoNum::Int32
      prec::Float32
      mwPerKm::Float32
      mvCl::Float32
end
system()=system(69,69.69,69.69,69.69)
####################################################################
mutable struct oya
   chichi::Tuple{Int32,Int32}
   haha::Tuple{Int32,Int32}
end
oya()=oya((0,0),(0,0))
####################################################################
mutable struct circuit
      binary::Array{Int8}
      decimal::Int32
      pcc::bus
      owpps::Array{bus}
      osss_owp::Array{bus}
      osss_mog::Array{bus}
      cost::Float32
      owp_MVcbls::Array{cbl}
      owp_HVcbls::Array{cbl}
      oss2oss_cbls::Array{cbl}
      pcc_cbls::Array{cbl}
      base_owp::bus
      oss_wind::wind
      oss_mva::Float32
      oyas::oya
      id::String
end
circuit()=circuit(Int8[],69,bus(),bus[],bus[],bus[],69.69,cbl[],cbl[],cbl[],cbl[],bus(),wind(),69.69,oya(),"id")
#########################################################################
mutable struct eqp_data
   cbls33kV::Array{Array{Float32, 1}, 1}
   cbls66kV::Array{Array{Float32, 1}, 1}
   cbls132kV::Array{Array{Float32, 1}, 1}
   cbls220kV::Array{Array{Float32, 1}, 1}
   cbls400kV::Array{Array{Float32, 1}, 1}
end
eqp_data()=eqp_data(Array{Array{Float32, 1}, 1}(),Array{Array{Float32, 1}, 1}(),Array{Array{Float32, 1}, 1}(),Array{Array{Float32, 1}, 1}(),Array{Array{Float32, 1}, 1}())
#########################################################################
mutable struct eez
      osss::Array{bus}
      owpps::Array{bus}
      pccs::Array{bus}
      buses::Int64
      sys::system
      finance::cstS_ks
      mv_circuits::Array{Array{circuit, 1}, 1}
      hv_circuits::Array{Array{circuit, 1}, 1}
      theta::Float32
      offset::Float32
      base::gps
      eqp_data::eqp_data
      mvc_pct::Float32
      hvc_pct::Float32
end
eez()=eez(bus[],bus[],bus[],0,system(),cstS_ks(),Array{Array{circuit, 1}, 1}(),Array{Array{circuit, 1}, 1}(),69.69,69.69,gps(),eqp_data(),0.0,0.0)
###################################################################
