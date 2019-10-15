mutable struct mutstruct1
    m2::mutstruct2
    mutstruct2() = (x = new(); x.m2 = x)
end

mutable struct mutstruct2
    a::Int64
    b::Int64
end
