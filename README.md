# SoftSurfaceDynamics.jl

Methods for integrating surface velocity field and for finding force equilibrium.

## What this package can do?

+ Create a sphere mesh subdivided from icosehadron. 
+ Calculate velocity due to a exerted force with `stokesvelocity` method. 
+ Calculate surface properties such as - normals, curvature, volume and area.
+ Stabilize tangential velocity components to keep mesh as good as possible uppon it evolution.

## What it can not do, but would be desirable to have in this package:

+ `meshpincher(msh)` 
+ `meshmerger(msh1,msh2)`
+ `patchmerger(msh1,msh2)`
+ `stabilize(msh)` a method for restructuring mesh by moving vertices in tangential plane and simplifying topology. (ElTopo.jl implementation in Julia)
+ `fixtopology(vertices,normals,topology)` a method for fixing connectivity and topology from meshes comming for example from marching cubes algorithm like in `Meshes.jl`. 

