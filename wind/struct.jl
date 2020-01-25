#this file describes the structure of wind objects

#wind object
mutable struct wind
      pu::Array{Float32}
      ce::Array{Float32}
      delta::Float32
      lf::Float32
end
wind()=wind(Float32[],Float32[],69.69,69.69)
