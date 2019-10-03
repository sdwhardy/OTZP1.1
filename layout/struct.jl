###################################################################
mutable struct line
      ymn::Float64
      ymax::Float64
      m::Float64
      b::Float64
end
line()=line(69.69,69.69,69.69,69.69)
###################################################################
mutable struct xy
      x::Float64
      y::Float64
end
xy()=xy(69.69,69.69)
###################################################################
mutable struct gps
      lat::Float64
      lng::Float64
end
gps()=gps(69.69,69.69)
###################################################################
mutable struct cntr
      gps::gps
      xy::xy
end
cntr()=cntr(gps(),xy())
###################################################################
mutable struct circle
      radius::Float64
      area::Float64
end
circle()=circle(69.69,69.69)
###################################################################
mutable struct bus
      mvas::Array{Float64}
      wnds::Array{String}
      inputs::Array{Int64}
      outputs::Array{Int64}
      cntr::cntr
      area::circle
      mv_area::circle
      num::Int64
      id::Int64
end
bus()=bus(Float64[],String[],Int64[],Int64[],cntr(),circle(),circle(),69,69)
###################################################################
mutable struct node
      cntr::cntr
end
node()=node(cntr())
###################################################################
mutable struct conductor
      head::bus
      tail::bus
      lngth::Float64
      mva::Float64#Only used for displaying solution
      num::Int64
      id::Int64
      #xfo_in
      #xfo_out
      #cbl
      #path[]
      #cost
      #kv
end
conductor()=conductor(bus(),bus(),69.69,69.69,69,69)
###################################################################
mutable struct edge
      head::node
      tail::node
      lngth::Float64
      build::Bool
end
edge()=edge(node(),node(),69.69,true)
####################################################################
mutable struct domain
      nodes::Array{node}
      edges::Array{edge}
      okNodes::Array{node}
      okEdges::Array{edge}
      nogoNodes::Array{node}
      nogoEdges::Array{edge}
end
domain()=domain(node[],edge[],node[],edge[],node[],edge[])
####################################################################
mutable struct nogo
      nodes::Array{node}
end
nogo()=nogo(node[])
####################################################################
mutable struct system
      nogoNum::Int64
      prec::Float64
end
system()=system(69,69.69)
####################################################################
mutable struct eez
      osss::Array{bus}
      owpps::Array{bus}
      pccs::Array{bus}
      buses::Array{bus}
      wp2pccs::Array{conductor}
      wp2osss::Array{conductor}
      oss2osss::Array{conductor}
      oss2pccs::Array{conductor}
      conductors::Array{conductor}
      bndryPnts::Array{node}
      nogos::Array{nogo}
      discretedom::domain
      sys::system
      theta::Float64
      offset::Float64
      base::gps
      id_count::Int64
end
eez()=eez(bus[],bus[],bus[],bus[],conductor[],conductor[],conductor[],conductor[],conductor[],node[],nogo[],domain(),system(),69.69,69.69,gps(),69)
###################################################################
