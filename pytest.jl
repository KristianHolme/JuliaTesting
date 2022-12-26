using Pkg
using PyCall
np = pyimport("numpy")
#Pkg.add("SciPy")
using SciPy
spsa = spatial
include("GaussQuad.jl")
using .GaussQuad
using .GetPlate
getPlate= GetPlate.getPlate
using PyPlot
p, tri, edge = getPlate(4)

