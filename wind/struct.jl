#this file describes the structure of wind objects

#wind object
mutable struct wind
      pu::Array{Float64}
      ce::Array{Float64}
      delta::Float64
      lf::Float64
end
wind()=wind(Float64[],Float64[],69.69,69.69)
