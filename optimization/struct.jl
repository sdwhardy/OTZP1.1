###################################################################
mutable struct ellipse
      y0::Float32
      x0::Float32
      ry::Float32
      rx::Float32
      alpha::Float32
end
ellipse()=ellipse(69.69,69.69,69.69,69.69,69.69)
###################################################################
mutable struct square
      ymx::Float32
      xmx::Float32
      ymn::Float32
      xmn::Float32
end
square()=square(69.69,69.69,69.69,69.69)
###################################################################
mutable struct constraints
      ellipses::Array{ellipse}
      #squares::Array{square}
end
#constraints()=constraints(ellipse[],square[])
constraints()=constraints(ellipse[])
###################################################################
