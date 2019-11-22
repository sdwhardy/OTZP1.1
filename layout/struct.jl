###################################################################
mutable struct line
      ymn::Float64
      ymx::Float64
      xmn::Float64
      xmx::Float64
      m_findx::Float64
      b_findx::Float64
      m_findy::Float64
      b_findy::Float64
end
line()=line(69.69,69.69,69.69,69.69,69.69,69.69,69.69,69.69)
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
      goal::Int64
      parent::node
      node()=(x=new();x.gps=gps();x.xy=xy();x.edges=Array{edge,1}();x.G_cost=Inf;x.H_cost=Inf;x.F_cost=Inf;x.num=69;x.openQ=false;x.closedQ=false; x.goal=0; x.parent=x)
end
#node()=(gps(),xy(),edge[],Inf,Inf,Inf,69,false,false,0,node())
###################################################################
mutable struct farm
      pos_height::Float64
      pos_width::Float64
      neg_width::Float64
      neg_height::Float64
      area::Float64
      nodes::Array{node}
      edges::Array{edge}
      ebnd::Array{line}
      wbnd::Array{line}
      sbnd::Array{line}
      nbnd::Array{line}
end
#farm()=farm(69.69,69.69,69.69,69.69,69.69,line(),line(),line(),line(),node[])
farm()=farm(69.69,69.69,69.69,69.69,69.69,node[],edge[],line[],line[],line[],line[])

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
###################################################################
#the structure used for a cable
mutable struct cbl
   mva::Float64
   length::Float64
   pth::Array{node}
   size::Float64
   cost::Float64
   num::Float64
   reliability::relia
   elec::elec
   costs::cbl_costs
end
cbl()=cbl(0.0,0.0,node[],0.0,0.0,0.0,relia(),elec(),cbl_costs())
###################################################################
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
###################################################################
mutable struct bus
      mva::Float64
      wnd::wind
      #inputs::Array{Int64}
      #outputs::Array{Int64}
      name::String
      node::node
      zone::farm
      mv_zone::farm
      kV::Int64
      num::Int64
      id::Int64
      xfmrs::Array{xfo}
end
bus()=bus(69.69,wind(),"sixty-nine",node(),farm(),farm(),69,69,69,xfo[])
###################################################################
#=mutable struct conductor
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
conductor()=conductor(bus(),bus(),69.69,69.69,69,69)=#
####################################################################
mutable struct domain
      nodes::Array{node}
      edges::Array{edge}
end
domain()=domain(node[],edge[])
####################################################################
mutable struct nogo
      nodes::Array{node}
      bndryPnts::Array{node}
      wbnd::Array{line}
      ebnd::Array{line}
      sbnd::Array{line}
      nbnd::Array{line}
end
nogo()=nogo(node[],node[],line[],line[],line[],line[])
####################################################################
mutable struct system
      nogoNum::Int64
      prec::Float64
      mwPerKm::Float64
end
system()=system(69,69.69,69.69)
####################################################################
mutable struct circuit
      binary::Array{Int8}
      decimal::Int64
      pcc::bus
      owpps::Array{bus}
      osss_owp::Array{bus}
      osss_mog::Array{bus}
      pths::Array{node}
      cost::Float64
      lengths::Array{Float64}
      owp_MVcbls::Array{cbl}
      owp_HVcbls::Array{cbl}
      oss2oss_cbls::Array{cbl}
      pcc_cbls::Array{cbl}
      parent_circ::Int64
      base_owp::bus
      oss_wind::wind
      oss_mva::Float64
      Qing::Bool
end
circuit()=circuit(Int8[],69,bus(),bus[],bus[],bus[],node[],69.69,Float64[],cbl[],cbl[],cbl[],cbl[],0,bus(),wind(),69.69,false)
#######################################################################################
mutable struct eez
      osss::Array{bus}
      owpps::Array{bus}
      pccs::Array{bus}
      buses::Array{bus}
      #wp2pccs::Array{conductor}
      #wp2osss::Array{conductor}
      #oss2osss::Array{conductor}
      #oss2pccs::Array{conductor}
      #conductors::Array{conductor}
      bndryPnts::Array{node}
      ebnd::Array{line}
      wbnd::Array{line}
      nbnd::Array{line}
      sbnd::Array{line}
      nogos::Array{nogo}
      discretedom::domain
      constrain::constraints
      sys::system
      finance::cstS_ks
      circuits::Array{circuit}
      theta::Float64
      offset::Float64
      base::gps
      id_count::Int64
end
#eez()=eez(bus[],bus[],bus[],bus[],conductor[],conductor[],conductor[],conductor[],conductor[],node[],nogo[],domain(),system(),cstS_ks(),69.69,69.69,gps(),69)
eez()=eez(bus[],bus[],bus[],bus[],node[],line[],line[],line[],line[],nogo[],domain(),constraints(),system(),cstS_ks(),circuit[],69.69,69.69,gps(),69)
###################################################################
