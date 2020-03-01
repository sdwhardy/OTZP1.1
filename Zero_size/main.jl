using DataFrames, XLSX, CSV, JLD2, FileIO
using StatsPlots, Plots
using Polynomials, TypedPolynomials, SpecialFunctions
using LinearAlgebra, Polyhedra, SetProg, MathOptInterface
using JuMP, Ipopt, Mosek, MosekTools

include("struct.jl")#
include("economics/functions.jl")#
include("layout/functions.jl")#



@time ocean=lof_layoutEez_basis()
@time ocean.owpps=lof_order2Pcc(ocean,ocean.pccs[2])
@time ocean=lof_layoutEez_expand(ocean,ocean.pccs[2])
