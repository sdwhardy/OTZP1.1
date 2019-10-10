#data structures used for equipment are specified in this file
mutable struct elec
   amp::Float64
   volt::Float64
   ohm::Float64
   farrad::Float64
   henry::Float64
   yc::Float64
   xl::Float64
end
elec()=elec(69.69,69.69,69.69,69.69,69.69,69.69,69.69)

#reliability structure
mutable struct relia
   fr::Float64
   mttr::Float64
   mc::Float64
end
relia()=relia(69.69,69.69,69.69)

#the structure for a transformer
mutable struct xfo
   mva::Float64
   num::Float64
   eta::Float64
   reliability::relia
   elec::elec
   costs::xfo_costs
end
xfo()=xfo(0.0,0.0,0.0,relia(),elec(),xfo_costs())

#the structure used for a cable
mutable struct cbl
   mva::Float64
   length::Float64
   size::Float64
   cost::Float64
   num::Float64
   reliability::relia
   elec::elec
   costs::cbl_costs
end
cbl()=cbl(0.0,0.0,0.0,0.0,0.0,relia(),elec(),cbl_costs())

#the structure used for a owpp (cable and xfm)
mutable struct owpp
   mva::Float64
   km::Float64
   cable::cbl
   xfm::xfo
   wp::wind
   costs::results
end
owpp()=owpp(0.0,0.0,cbl(),xfo(),wind(),results())
