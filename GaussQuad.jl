module GaussQuad
    using LinearAlgebra

    function refToab(y, a, b)
        return ((b-a)*y + b + a)/2
    end

    function Area(p)::Float64
        p1, p2, p3 = p[1], p[2], p[3]
        return 0.5*norm( (p2[1]-p1[1])*(p3[2]-p1[2]) - (p2[2]-p1[2])*(p3[1]-p1[1]) )
    end

    function BaryToCart(p1, p2, p3, b)
        S = p1*b[1] + p2*b[2] + p3*b[3]
        return S
    end

    function quadrature1D(a, b, Nq, g)
        nodesN1 = [0.0]
        nodesN2 = [-1/sqrt(3), 1/sqrt(3) ]
        nodesN3 = [-sqrt(3/5), 0, sqrt(3/5)]
        nodesN4 = [-sqrt((3+2*sqrt(6/5))/7), -sqrt((3-2*sqrt(6/5))/7),
                            sqrt((3-2*sqrt(6/5))/7), sqrt((3+2*sqrt(6/5))/7)]
        nodesN5 = [-1/3 * sqrt(5+2*sqrt(10/7)), -1/3 * sqrt(5-2*sqrt(10/7)), 
                            0, 1/3 * sqrt(5-2*sqrt(10/7)), 1/3 * sqrt(5+2*sqrt(10/7))]

        weightsN1 = [2]
        weightsN2 = [1, 1]
        weightsN3 = [5/9, 8/9, 5/9]
        weightsN4 = [(18 - sqrt(30))/36, (18 + sqrt(30))/36, 
                            (18+sqrt(30))/36, (18-sqrt(30))/36]
        weightsN5 = [(322-13*sqrt(70))/900, (322+13*sqrt(70))/900, 128/225, 
                            (322+13*sqrt(70))/900, (322-13*sqrt(70))/900]

        nodes = [nodesN1, nodesN2, nodesN3, nodesN4, nodesN5]
        weights = [weightsN1, weightsN2, weightsN3, weightsN4, weightsN5]

        chosenWeights = weights[Nq]
        chosenNodes = nodes[Nq]


        I = (b-a)/2 * chosenWeights'*(g.(refToab.(chosenNodes, a, b)))

        return I
    end


    function quadrature2D(p, Nq, g)
        p1, p2, p3 = p[1], p[2], p[3]
        
        
        nodesN1 = [[1/3, 1/3, 1/3]]
        nodesN3 = [[1/2, 1/2, 0], [1/2, 0, 1/2], [0, 1/2, 1/2]]
        nodesN4 = [[1/3, 1/3, 1/3], [3/5, 1/5, 1/5], [1/5, 3/5, 1/5], [1/5, 1/5, 3/5]]
        
        weightsN1 = [1]
        weightsN3 = [1/3, 1/3, 1/3]
        weightsN4 = [-9/16, 25/48, 25/48, 25/48]
        
        nodes = [nodesN1,NaN, nodesN3, nodesN4]
        weights = [weightsN1, NaN, weightsN3, weightsN4]
        
        i = Nq
        chosenNodes, chosenWeights = nodes[i], weights[i]
        
        gValues = g.( map(x -> BaryToCart(p1, p2, p3, x), chosenNodes ) )
        A = Area(p)
        
        I = A* chosenWeights'*gValues
    end


    function LineQuadrature(a, b, Nq, g)
        length = norm(b-a)
        
        f(x) = g( a*(1-x) + x*b )
        
        I = length*quadrature1D(0, 1, Nq, f)
        return I
    end
end