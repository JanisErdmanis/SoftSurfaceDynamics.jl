### Some methods for velocity calculation for a system driven by a surface force 
"""
Interface velocity in a viscous liquid where regions seperated by interface do have the same viscousity. If γ is passed the curvattureless algorithm is used to take into account force due to surface tension. Returns velocity field projected on vertex normals. Usefull for equilibrium calculations.
"""
function stokesvelocity(points,normals,faces,forcen,etaP,gammap)
    vareas = vertexareas(points,faces)
    velocityn = zeros(Float64,length(points))

    for xkey in 1:length(points)

        x = points[xkey]
        nx = normals[xkey]
        fx = forcen[xkey]
                
        s = 0
        for ykey in 1:length(points)
            if ykey==xkey
                continue
            end

            y = points[ykey]
            ny = normals[ykey]
            fy = forcen[ykey]

            ### I will need to check a missing 2
            s += vareas[ykey]*1 ./8/pi/etaP* dot(y-x,nx+ny)/norm(y-x)^3*(1-3*dot(y-x,nx)*dot(y-x,ny)/norm(y-x)^2) * gammap

            ### ????????
            s += vareas[ykey]*1 ./8/pi/etaP* ( dot(nx,ny)/norm(x-y) + dot(nx,x -y)*dot(ny,x-y)/norm(x-y)^3 )*(fy - fx)
        end
        velocityn[xkey] = s
    end
    return velocityn
end


function stokesvelocity(points,normals,faces,forcen,etaP)

    vareas = vertexareas(points,faces)
    velocityn = zeros(Float64,size(points,2))
    
    for xkey in 1:length(points)

        x = points[xkey]
        nx = normals[xkey]
        fx = forcen[xkey]
                
        s = 0
        for ykey in 1:length(points)
            if ykey==xkey
                continue
            end

            y = points[ykey]
            ny = normals[ykey]
            fy = forcen[ykey]
            
            s += vareas[ykey]*1 ./8/pi/etaP* ( dot(nx,ny)/norm(x-y) + dot(nx,x -y)*dot(ny,x-y)/norm(x-y)^3 )*(fy - fx)
        end

        velocityn[xkey] = s
    end

    return velocityn
end
