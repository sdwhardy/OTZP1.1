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
mutable struct edge
      head::Int64
      tail::Int64
      lngth::Float64
end
edge()=edge(69,69,69.69)
###################################################################
mutable struct node
      gps::gps
      xy::xy
      edges::Array{edge}
      G_cost::Float64
      H_cost::Float64
      F_cost::Float64
      num::Int64
      openQ::Bool
      closedQ::Bool
      parent::node
      node()=(x=new();x.gps=gps();x.xy=xy();x.edges=Array{edge,1}();x.G_cost=Inf;x.H_cost=Inf;x.F_cost=Inf;x.num=69;x.openQ=false;x.closedQ=false; x.parent=x)
end
###################################################################
mutable struct farm
      pos_height::Float64
      pos_width::Float64
      neg_width::Float64
      neg_height::Float64
      area::Float64
      pnts::Array{node}
      edges::Array{edge}
end
#farm()=farm(69.69,69.69,69.69,69.69,69.69,line(),line(),line(),line(),node[])
farm()=farm(69.69,69.69,69.69,69.69,69.69,node[],edge[])
###################################################################
mutable struct bus
      mvas::Array{Float64}
      wnds::Array{wind}
      inputs::Array{Int64}
      outputs::Array{Int64}
      name::String
      node::node
      zone::farm
      mv_zone::farm
      num::Int64
      id::Int64
end
bus()=bus(Float64[],wind[],Int64[],Int64[],"sixty-nine",node(),farm(),farm(),69,69)
###################################################################
mutable struct conductor
      head::bus
      tail::bus
      lngth::Float64
      mva::Float64
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
####################################################################
mutable struct domain
      nodes::Array{node}
      edges::Array{edge}
end
domain()=domain(node[],edge[])
####################################################################
mutable struct nogo
      nodes::Array{node}
end
nogo()=nogo(node[])
####################################################################
mutable struct system
      nogoNum::Int64
      prec::Float64
      mwPerKm::Float64
end
system()=system(69,69.69,69.69)
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
      finance::cstS_ks
      theta::Float64
      offset::Float64
      base::gps
      id_count::Int64
end
eez()=eez(bus[],bus[],bus[],bus[],conductor[],conductor[],conductor[],conductor[],conductor[],node[],nogo[],domain(),system(),cstS_ks(),69.69,69.69,gps(),69)
###################################################################
