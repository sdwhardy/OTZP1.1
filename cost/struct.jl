#=
This file defines the structure of any objects associated with costs
=#

#the structure of costs for a cable
mutable struct cbl_costs
   qc::Float64
   cbc::Float64
   rlc::Float64
   cm::Float64
   eens::Float64
   ttl::Float64
end
cbl_costs()=cbl_costs(0.0,0.0,0.0,0.0,0.0,0.0)

#the structure of costs for a transformers
mutable struct xfo_costs
   cpx::Float64
   tlc::Float64
   cm::Float64
   eens::Float64
   ttl::Float64
end
xfo_costs()=xfo_costs(0.0,0.0,0.0,0.0,0.0)

#cost components and totals calculated for a OWPP object
mutable struct results
     cpx::Float64
     loss::Float64
     opex::Float64
     ttl::Float64
end
results()=results(0.0,0.0,0.0,0.0)

#an object that contains all cost factors used in the calculations
mutable struct cstS_ks
   FC_ac::Float64
   FC_dc::Float64
   dc::Float64
   f_ct::Float64
   p_ct::Float64
   c_ct::Float64
   Qc_oss::Float64
   Qc_pcc::Float64
   life::Float64
   T_op::Float64
   E_op::Float64
   cf::Float64
end
cstS_ks()=cstS_ks(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
