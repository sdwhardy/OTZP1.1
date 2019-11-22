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
