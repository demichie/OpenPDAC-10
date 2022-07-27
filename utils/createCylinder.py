import meshzoo
import numpy as np
from stl import mesh

xbase_center = 0.0
ybase_center = 0.0
zbase_center = 800.0

radius = 10.0
height = 1200.0

n_subdivide = 50
vertices, faces = meshzoo.tetra_sphere(n_subdivide)

vertices, faces = meshzoo.tube(length=height, radius=radius, n=50)

vertices = vertices + [xbase_center,ybase_center,0.5*height+zbase_center]

# Create the mesh
cylinder = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        cylinder.vectors[i][j] = vertices[f[j],:]

# Write the mesh to file "cube.stl"
cylinder.save('../constant/triSurface/cylinder.stl')

