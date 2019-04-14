using Optim
using Parameters
using LinearAlgebra
using GeometryTypes

### If C=0 we have original Zinchenko 1997 stabilization scheme
@with_kw struct KineticStabilizer
    C::AbstractFloat = 0.4
    ftol::AbstractFloat = 1e-6
end

function F(v,points,faces,zc::KineticStabilizer)

    Cp = zc.C
    
    s = 0    
    for ti in 1:length(faces)
        v1,v2,v3 = faces[ti]

        a = norm(points[v2] - points[v1])
        b = norm(points[v3] - points[v2])
        c = norm(points[v1] - points[v3])

        xava = dot(points[v2]-points[v1],v[v2]-v[v1])
        xbvb = dot(points[v3]-points[v2],v[v3]-v[v2])
        xcvc = dot(points[v1]-points[v3],v[v1]-v[v3])

        ### This part is responsible for first part in the sum
        s += v2>v1 ? 4 * xava^2 : 0
        s += v3>v2 ? 4 * xbvb^2 : 0
        s += v1>v3 ? 4 * xcvc^2 : 0

        ### This part is needed for ordinary things
        Cdelta = 1/4*sqrt(1 - 2* (a^4 + b^4 + c^4)/(a^2 + b^2 + c^2)^2)
        A = 1/4/Cdelta/(a^2 + b^2 + c^2)^3 * ( a^2*(a^2 + b^2 + c^2) - a^4 - b^4 - c^4)
        B = 1/4/Cdelta/(a^2 + b^2 + c^2)^3 * ( b^2*(a^2 + b^2 + c^2) - a^4 - b^4 - c^4)
        C = 1/4/Cdelta/(a^2 + b^2 + c^2)^3 * ( c^2*(a^2 + b^2 + c^2) - a^4 - b^4 - c^4)
        DCdelta = -A*xava - B*xbvb - C*xcvc

        s += Cp*DCdelta^2/Cdelta^2
    end

    return s
end

function gradF!(v,points,faces,storage,zc::KineticStabilizer)

    Cp = zc.C

    for i in 1:length(storage)
        storage[i] = Point(0,0,0)
    end
    
    for i in 1:length(faces)
        v1,v2,v3 = faces[i]
        a = norm(points[v2] - points[v1])
        b = norm(points[v3] - points[v2])
        c = norm(points[v1] - points[v3])

        Cdelta = 1/4*sqrt(1 - 2 * (a^4 + b^4 + c^4)/(a^2 + b^2 + c^2)^2)
        A = 1/4/Cdelta/(a^2 + b^2 + c^2)^3 * ( a^2*(a^2 + b^2 + c^2) - a^4 - b^4 - c^4)
        B = 1/4/Cdelta/(a^2 + b^2 + c^2)^3 * ( b^2*(a^2 + b^2 + c^2) - a^4 - b^4 - c^4)
        C = 1/4/Cdelta/(a^2 + b^2 + c^2)^3 * ( c^2*(a^2 + b^2 + c^2) - a^4 - b^4 - c^4)
        
        xa = points[v2] - points[v1]
        xb = points[v3] - points[v2]
        xc = points[v1] - points[v3]

        xava = dot(points[v2]-points[v1],v[v2]-v[v1])
        xbvb = dot(points[v3]-points[v2],v[v3]-v[v2])
        xcvc = dot(points[v1]-points[v3],v[v1]-v[v3])

        DCdelta = -A*xava - B*xbvb - C*xcvc

        A_ = 2*Cp*DCdelta*A/Cdelta^2 
        B_ = 2*Cp*DCdelta*B/Cdelta^2 
        C_ = 2*Cp*DCdelta*C/Cdelta^2 
        
        storage[v1] += (A_ - 8*xava)*xa - C_*xc
        storage[v2] += (B_ - 8*xbvb)*xb - A_*xa
        storage[v3] += (C_ - 8*xcvc)*xc - B_*xb
    end
end

function stabilize!(v,points,faces,n,zc::KineticStabilizer)
    #h2 = hvec(points,faces,n)

    ### Let's try to make it simpler 
    
    function f(x::Vector)
        # vv = reshape(x,3,div(length(x),3))
        return F(x,points,faces,zc)
    end

    function g!(storage::Vector,x::Vector)

        # vv = reshape(x,size(points)...)
        # gradf = reshape(storage,size(points)...)
        gradF!(x,points,faces,storage,zc)
        for i in 1:size(v,2)
            P = I - n[i]*n[i]'
            storage[i] = P*storage[i]
        end
    end
    
    res = optimize(f,g!,v,method=Optim.ConjugateGradient(),ftol=zc.ftol)

    ### One needs to be a little bit carefull about what types would Julia return
    return res
    #v[:,:] = reshape(res.minimum, size(v)...)[:,:]
end
