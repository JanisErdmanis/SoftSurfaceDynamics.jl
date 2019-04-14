import SoftSurfaceDynamics

msh = SoftSurfaceDynamics.unitsphere(3)

normals = SoftSurfaceDynamics.normals(msh.vertices,msh.faces)
curvatures = SoftSurfaceDynamics.meancurvatures(msh.vertices,msh.faces)
vareas = SoftSurfaceDynamics.vertexareas(msh.vertices,msh.faces)

volume = SoftSurfaceDynamics.surfacevolume(msh.vertices,msh.faces)
