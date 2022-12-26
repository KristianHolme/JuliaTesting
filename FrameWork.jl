include("./GaussQuad.jl")
using.GaussQuad
include("./GetPlate.jl")
using .GetPlate
getPlate= GetPlate.getPlate
using PyCall
using Conda
#Conda.add("matplotlib")

mpl = pyimport("matplotlib")


function C(k, p, tri)
    """
    returns a 3x3 matrix where C, where each column contains
    c^k_alpha, c^k_{x alpha}, and c^k_{y alpha}
    
    parameters:
        k: element number k, k>=0
        p: points in triangulation, from getdisc
        tri: triangulation, from getdisc
        
    """
    p1 = p[tri[k, 1]]
    p2 = p[tri[k, 2]]
    p3 = p[tri[k, 3]]
    A = [[ 1, p1[0], p1[1] ], 
        [ 1, p2[0], p2[1] ], 
        [ 1, p3[0], p3[1] ] ]
    C = inv(A)
    return C
end

m = [1 0 0 
     0 1 0
     0 0 1]
     
p, tri, edge = getPlate(4)

