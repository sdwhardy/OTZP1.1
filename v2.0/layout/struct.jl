mutable struct system
      nogoNum::Int32
      prec::Float32
      mwPerKm::Float32
      mvCl::Float32
end
system()=system(69,69.69,69.69,69.69)

mutable struct node
      gps::gps
      xy::xy
      num::Int32
end
node()=node(gps(),xy(),69)

mutable struct bus
      mva::Float32
      wnd::wind
      name::String
      node::node
      mv_zone::Float32
      kv2pcc::Float32
      kV::Float32
      xfmrs::Array{transformer}
      conv::Array{converter}
      plat::Array{platform}
end
bus()=bus(69.69,wind(),"sixty-nine",node(),69.69,69.69,69.69,transformer[],converter[],platform[])

mutable struct circuit
      binary::Array{Int8}
      decimal::Int32
      pcc::bus
      owpps::Array{bus}
      oss::Array{bus}
      mog::Array{bus}
      cost::Float32
      MVcbls::Array{cable}
      HVcbls::Array{cable}
      O2Ocbls::Array{cable}
      PCCcbls::Array{cable}
      base::bus
      wnd::wind
      mva::Float32
      id::String
end
circuit()=circuit(Int8[],69,bus(),bus[],bus[],bus[],69.69,cable[],cable[],cable[],cable[],bus(),wind(),69.69,"id")

mutable struct eez
      owpps::Array{bus}
      pcc::bus
      sys::system
      num::Int32#total number of nodes in the system
      mv_circuits::Array{Array{circuit, 1}, 1}
      hv_circuits::Array{Array{circuit, 1}, 1}
      offset::xy
      database::Dict{String, Dict{String,Any}}
end
eez()=eez(bus[],bus(),system(),0,Array{Array{circuit, 1}, 1}(),Array{Array{circuit, 1}, 1}(),xy(),Dict{String, Dict{String,Any}}())
