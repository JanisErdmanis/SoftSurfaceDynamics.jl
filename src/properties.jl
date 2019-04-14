using GeometryTypes
using SurfaceTopology
using LinearAlgebra

function quadraticform(vects,vnormal)
    
    Lx = [0 0 0; 0 0 -1; 0 1 0]
    Ly = [0 0 1; 0 0 0; -1 0 0]
    Lz = [0 -1 0; 1 0 0; 0 0 0]

    d = [0,0,1] + vnormal
    d /= norm(d)
        
    Ln = d[1]*Lx + d[2]*Ly + d[3]*Lz
    R = exp(pi*Ln)

    vects = copy(vects)
    for vj in 1:length(vects)
        vects[vj] = R*vects[vj]
    end

    ### Construction of the system
    A = Array{Float64}(undef,3,3)
    B = Array{Float64}(undef,3)

    vects_norm2 = Array{Float64}(undef,length(vects))
    for vj in 1:length(vects)
       vects_norm2[vj] = norm(vects[vj])^2
    end

    A[1,1] = sum((v[1]^4 for v in vects) ./ vects_norm2)
    A[1,2] = sum((v[1]^3*v[2] for v in vects) ./ vects_norm2)
    A[1,3] = sum((v[1]^2*v[2]^2 for v in vects) ./ vects_norm2)
    A[2,1] = A[1,2]
    A[2,2] = A[1,3]
    A[2,3] = sum( (v[2]^3*v[1] for v in vects) ./vects_norm2)
    A[3,1] = A[1,3]
    A[3,2] = A[2,3]
    A[3,3] = sum((v[2]^4 for v in vects) ./vects_norm2)

    
    B[1] = sum((v[3]*v[1]^2 for v in vects) ./vects_norm2)
    B[2] = sum((v[1]*v[2]*v[3] for v in vects) ./vects_norm2)
    B[3] = sum((v[2]^2*v[3] for v in vects) ./vects_norm2)

    C,D,E = A\B
    return C,D,E
end

function meancurvature(points,topology)
    curvatures = Array{Float64}(undef,length(points))
    for v in 1:length(points)

        s = Point(0,0,0)
        for (v1,v2) in EdgeRing(v,topology)
            s += cross(points[v2],points[v1])
        end
        normal = s ./ norm(s)

        vring = collect(VertexRing(v,topology))
        vects = [points[vi] - points[v] for vi in vring]

        C,D,E = quadraticform(vects,normal)

        A = [C D/2;D/2 E]
        k1,k2 = eigvals(-A)
        H = (k1 + k2)/2

        ### Multiplier by is 2 an empirical fix
        curvatures[v] = 2*H 
    end
    return curvatures
end

function normals(vertices,topology)
    n = Point{3,Float64}[]
    for v in 1:length(vertices)
        s = Point(0,0,0)
        for (v1,v2) in EdgeRing(v,topology)
            s += cross(vertices[v2],vertices[v1])
        end
        normal = s ./ norm(s)
        push!(n,normal)
    end
    return n
end

function vertexareas(points,topology)
    vareas = zeros(Float64,length(points))
    for face in Faces(topology)
        v1,v2,v3 = face
        area = norm(cross(points[v2]-points[v1],points[v3]-points[v1])) /2
        vareas[v1] += area/3
        vareas[v2] += area/3
        vareas[v3] += area/3
    end
    return vareas
end

function surfacevolume(points,topology)
    # Calculate face normal
    # Calculate projected area
    # Calculate ordinary volume 
    # Calculate volume between projected and real area
    # (+) if normal is outwards
    
    normal0 = [0,0,1]
    s = 0
    for face in Faces(topology)
        y1 = points[face[1]]
        y2 = points[face[2]]
        y3 = points[face[3]]

        normaly = cross(y2-y1,y3-y1)
        normaly /= norm(normaly)

        area = norm(cross(y2-y1,y3-y1))/2
        areaproj = dot(normaly,normal0)*area
        volume = dot(y1 + y2 + y3,normal0)/3*areaproj

        s += volume
    end

    return s
end



