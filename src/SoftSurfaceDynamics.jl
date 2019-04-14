module SoftSurfaceDynamics
### Tools necessary to compute and integrate surface evolution

### Things which could belong to this package

# meshpincher(msh)
# meshmerger(msh1,msh2) 
# patchmerger(msh1,msh2)
# curvature calculation [Done]
# fixtopology(msh)
# meshstabilizer

include("properties.jl")
include("velocity.jl")
include("meshes.jl")
include("kineticstabilizer.jl")

export unitsphere, meancurvatures, normals, surfacevolume, vertexareas, stokesvelocity
# export KineticStabilizer, stabilize!

end # module
