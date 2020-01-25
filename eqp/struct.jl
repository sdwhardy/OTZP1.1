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

#reliability structure
mutable struct relia
   fr::Float32
   mttr::Float32
   mc::Float32
end
relia()=relia(69.69,69.69,69.69)
