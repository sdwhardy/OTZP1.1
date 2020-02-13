###################################################################
mutable struct line
      ymn::Float32
      ymx::Float32
      xmn::Float32
      xmx::Float32
      m_findx::Float32
      b_findx::Float32
      m_findy::Float32
      b_findy::Float32
end
line()=line(69.69,69.69,69.69,69.69,69.69,69.69,69.69,69.69)
###################################################################
mutable struct xy
      x::Float32
      y::Float32
end
xy()=xy(69.69,69.69)
###################################################################
mutable struct gps
      lat::Float32
      lng::Float32
end
gps()=gps(69.69,69.69)
###################################################################
mutable struct edge
      head::Int32
      tail::Int32
      lngth::Float32
end
edge()=edge(69,69,69.69)
###################################################################
mutable struct node
      gps::gps
      xy::xy
      edges::Array{edge}
      G_cost::Float32
      H_cost::Float32
      F_cost::Float32
      num::Int32
      openQ::Bool
      closedQ::Bool
      goal::Int32
      parent::node
      node()=(x=new();x.gps=gps();x.xy=xy();x.edges=Array{edge,1}();x.G_cost=Inf;x.H_cost=Inf;x.F_cost=Inf;x.num=69;x.openQ=false;x.closedQ=false; x.goal=0; x.parent=x)
end
#node()=(gps(),xy(),edge[],Inf,Inf,Inf,69,false,false,0,node())
###################################################################
mutable struct farm
      pos_height::Float32
      pos_width::Float32
      neg_width::Float32
      neg_height::Float32
      area::Float32
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
   mva::Float32
   num::Float32
   hv::Float32
   lv::Float32
   eta::Float32
   reliability::relia
   elec::elec
   costs::xfo_costs
end
xfo()=xfo(0.0,0.0,0.0,0.0,0.0,relia(),elec(),xfo_costs())
###################################################################
#the structure used for a cable
mutable struct cbl
   mva::Float32
   length::Float32
   pth::Array{node}
   size::Float32
   cost::Float32
   num::Float32
   reliability::relia
   elec::elec
   costs::cbl_costs
end
cbl()=cbl(0.0,0.0,node[],0.0,0.0,0.0,relia(),elec(),cbl_costs())
###################################################################
mutable struct bus
      mva::Float32
      wnd::wind
      #inputs::Array{Int32}
      #outputs::Array{Int32}
      name::String
      node::node
      zone::farm
      mv_zone::farm
      kV::Int32
      num::Int32
      id::Int32
      xfmrs::Array{xfo}
      base_cost::Float32
end
bus()=bus(69.69,wind(),"sixty-nine",node(),farm(),farm(),69,69,69,xfo[],0)
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
      nogoNum::Int32
      prec::Float32
      mwPerKm::Float32
      mvCl::Float32
end
system()=system(69,69.69,69.69,69.69)
####################################################################
mutable struct circuit
      binary::Array{Int8}
      decimal::Int32
      pcc::bus
      owpps::Array{bus}
      osss_owp::Array{bus}
      osss_mog::Array{bus}
      pths::Array{node}
      cost::Float32
      lengths::Array{Float32}
      owp_MVcbls::Array{cbl}
      owp_HVcbls::Array{cbl}
      oss2oss_cbls::Array{cbl}
      pcc_cbls::Array{cbl}
      parent_circ::Int32
      base_owp::bus
      oss_wind::wind
      oss_mva::Float32
      Qing::Bool
end
circuit()=circuit(Int8[],69,bus(),bus[],bus[],bus[],node[],69.69,Float32[],cbl[],cbl[],cbl[],cbl[],0,bus(),wind(),69.69,false)
#######################################################################################
mutable struct eez
      osss::Array{bus}
      owpps::Array{bus}
      pccs::Array{bus}
      buses::Array{bus}
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
      theta::Float32
      offset::Float32
      base::gps
      id_count::Int32
      yaxisMajor::Bool
end
eez()=eez(bus[],bus[],bus[],bus[],node[],line[],line[],line[],line[],nogo[],domain(),constraints(),system(),cstS_ks(),circuit[],69.69,69.69,gps(),69,false)
###################################################################


########################################## depricated ##########################

###################################################################
#=mutable struct conductor
      head::bus
      tail::bus
      lngth::Float32
      mva::Float32
      num::Int32
      id::Int32
      #xfo_in
      #xfo_out
      #cbl
      #path[]
      #cost
      #kv
end
conductor()=conductor(bus(),bus(),69.69,69.69,69,69)=#

###################################################################
#the structure used for a owpp (cable and xfm)
#=mutable struct owpp
   mva::Float32
   km::Float32
   cable::cbl
   xfm::xfo
   wp::wind
   costs::results
end
owpp()=owpp(0.0,0.0,cbl(),xfo(),wind(),results())=#
