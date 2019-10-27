###################################################################
mutable struct ellipse
      y0::Float64
      x0::Float64
      ry::Float64
      rx::Float64
end
ellipse()=ellipse(69.69,69.69,69.69,69.69)
###################################################################
mutable struct square
      ymx::Float64
      xmx::Float64
      ymn::Float64
      xmn::Float64
end
square()=square(69.69,69.69,69.69,69.69)
###################################################################
mutable struct constraints
      ellipses::Array{ellipse}
      squares::Array{square}
end
constraints()=constraints(ellipse[],square[])
###################################################################
