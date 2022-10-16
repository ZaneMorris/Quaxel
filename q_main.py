#%% Intersection Testing
from Q_subrout1 import *
from Q_subrout2 import *

path = r'C:\Users\MDO-snow\Documents\05_SMSDProjects\05_Thesis\01_Code\01_QUAXEL\TestGeometries\STL'
filename = 'MasterJaw.stl'
facets, bounds = read_stl(path, filename)

# Create qtree with bounds
v_mesh = qtree(facets,th_nd=2, th_br=4,boundary=bounds)

# Generate quadtree stucture according to the conditions and stl file
v_mesh.subdivide()
print(len(v_mesh.root.tris))

# Draw Quadtree
#v_mesh.draw_quad()

# Create the enclosing voxel mesh
el = v_mesh.voxel_mesh(40)

# Determine which voxels are filled, and which are void
v_mesh.voxelize(el)
print(v_mesh.vox_arr[:,5,7])

# Plot the resulting voxel mesh
print('now plotting')
#v_mesh.draw_temp(el, alpha=0.9, color='b')
v_mesh.draw2()