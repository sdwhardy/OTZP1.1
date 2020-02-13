#=
This file defines the structure of any objects associated with costs
=#

#the structure of costs for a cable
mutable struct cbl_costs
   qc::Float32
   cbc::Float32
   rlc::Float32
   cm::Float32
   eens::Float32
   ttl::Float32
end
cbl_costs()=cbl_costs(0.0,0.0,0.0,0.0,0.0,0.0)

#the structure of costs for a transformers
mutable struct xfo_costs
   cpx::Float32
   tlc::Float32
   cm::Float32
   eens::Float32
   ttl::Float32
end
xfo_costs()=xfo_costs(0.0,0.0,0.0,0.0,0.0)

#cost components and totals calculated for a OWPP object
mutable struct results
     cpx::Float32
     loss::Float32
     opex::Float32
     ttl::Float32
end
results()=results(0.0,0.0,0.0,0.0)

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
end
cstS_ks()=cstS_ks(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
