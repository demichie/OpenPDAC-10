import meshzoo
import numpy as np
from stl import mesh

radius = 100.0
x_center = 0.0
y_center = 0.0
z_center = 1100.0


n_subdivide = 50
vertices, faces = meshzoo.tetra_sphere(n_subdivide)

vertices = vertices * radius

vertices = vertices + [x_center,y_center,z_center]

# Create the mesh
sphere = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
for i, f in enumerate(faces):
    for j in range(3):
        sphere.vectors[i][j] = vertices[f[j],:]

# Write the mesh to file "cube.stl"
sphere.save('../constant/triSurface/sphere.stl')

