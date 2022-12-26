module GetPlate

    using Pkg
    using PyCall
    np = pyimport("numpy")
    #Pkg.add("SciPy")
    using SciPy
    spsa = spatial


    function getPlate(N)
        # Defining auxiliary variables.
        L = np.linspace(-1,1,N)
        Y,X = np.meshgrid(L,L)
        x = np.ravel(np.transpose(X))
        y = np.ravel(np.transpose(Y))
        # Generating nodal points.
        p = np.zeros((N^2,2))
        for i in 1:N^2
            p[i,1] = x[i]
            p[i,2] = y[i]
        end
        # Generating elements.
        mesh = spsa.Delaunay(p)
        tri = mesh.simplices .+ 1

        # Generating nodal points on outer edge.
        south = np.array([np.arange(1,N),np.arange(2,N+1)])
        east = np.array([np.arange(N,N^2-N+1,N),np.arange(2*N,N^2+1,N)])
        north = np.array([np.arange(N^2,N^2-N+1,-1),np.arange(N^2-1,N^2-N,-1)])
        west = np.array([np.arange(N^2-N+1,N-1,-N),np.arange(N^2-2*N+1,0,-N)])
        L1 = np.shape(south)[1]
        L2 = np.shape(east)[1]
        L3 = np.shape(west)[1]
        L4 = np.shape(north)[1]
        edge = zeros(Int32, L1+L2+L3+L4,2)
        for i in 1:L1
            edge[i,1] = south[1,i]
            edge[i,2] = south[2,i]
        end
        for i in (L1+1):(L1+L2)
            edge[i,1] = east[1,i-L1]
            edge[i,2] = east[2,i-L1]
        end
        for i in (L1+L2+1):(L1+L2+L3)
            edge[i,1] = north[1,i-L1-L2]
            edge[i,2] = north[2,i-L1-L2]
        end
        for i in (L1+L2+L3+1):(L1+L2+L3+L4)
            edge[i,1] = west[1,i-L1-L2-L3]
            edge[i,2] = west[2,i-L1-L2-L3]
        end
        return p,tri,edge
    end

end
