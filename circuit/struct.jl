mutable struct circuit
      binary::Array{Int8}
      decimal::Int64
      pcc::bus
      pcc_pth::node
      pcc_cst::Float64
      pcc_length::Float64
      owpps::Array{bus}
      owp_pths::Array{node}
      owp_csts::Array{Float64}
      owp_lengths::Array{Float64}
      base_owp::bus
end
circuit()=circuit(Int8[],Int64,bus(),node(),69.69,69.69,bus[],node[],Float64[],Float64[],bus())
