from math import sqrt
import numpy as np
from pymicro.view.scene3d import Scene3D
from pymicro.view.vtk_utils import *
from pymicro.crystal.lattice import HklPlane
from pymicro.crystal.microstructure import Orientation

'''
Create a 3d scene with a hexagonal crystal lattice.
Hkl planes are added to the lattice and displayed.
'''
# Create the 3D scene
base_name = os.path.splitext(__file__)[0]
s3d = Scene3D(display=True, ren_size=(800, 800), name=base_name)

# hexagonal lattice for Ti7Al
a = 0.2931  # nm
c = 0.4694  # nm
l = Lattice.hexagonal(a, c)

"""
  0   0   0   1  basal 1
"""
print('basal plane:')
print(HklPlane.four_to_three_indices(0, 0, 0, 1))  # basal
basal = HklPlane(0, 0, 1, lattice=l)
hkl_planes = [basal]
"""
  0   1  -1   0  prism 1
 -1   1   0   0  prism 2
 -1   0   1   0  prism 3
"""
print('prismatic planes:')
print(HklPlane.four_to_three_indices(0, 1, -1, 0))  # prismatic 1
print(HklPlane.four_to_three_indices(-1, 1, 0, 0))  # prismatic 2
print(HklPlane.four_to_three_indices(-1, 0, 1, 0))  # prismatic 3
prism_1 = HklPlane(0, 1, 0, lattice=l)
prism_2 = HklPlane(-1, 1, 0, lattice=l)
prism_3 = HklPlane(-1, 0, 0, lattice=l)
#hkl_planes = [prism_1, prism_2, prism_3]
hkl_planes.extend([prism_1, prism_2, prism_3])
"""
  0   1  -1   1  pyr 1
 -1   1   0   1  pyr 2
 -1   0   1  -1  pyr 3
 -1   0   1   1  pyr 4
 -1   1   0  -1  pyr 5
  0  -1   1   1  pyr 6
"""
print('pyramidal planes:')
print(HklPlane.four_to_three_indices(0, 1, -1, 1))  # prismatic 1
print(HklPlane.four_to_three_indices(-1, 1, 0, 1))  # prismatic 2
print(HklPlane.four_to_three_indices(-1, 0, 1, -1))  # prismatic 3
print(HklPlane.four_to_three_indices(-1, 0, 1, 1))  # prismatic 1
print(HklPlane.four_to_three_indices(-1, 1, 0, -1))  # prismatic 2
print(HklPlane.four_to_three_indices(0, -1, 1, 1))  # prismatic 3
pyr_1 = HklPlane(0, 1, 1, lattice=l)
pyr_2 = HklPlane(-1, 1, 1, lattice=l)
pyr_3 = HklPlane(-1, 0, -1, lattice=l)
pyr_4 = HklPlane(-1, 0, 1, lattice=l)
pyr_5 = HklPlane(-1, 1, -1, lattice=l)
pyr_6 = HklPlane(0, -1, 1, lattice=l)
#hkl_planes = [pyr_1, pyr_2, pyr_3, pyr_4, pyr_5, pyr_6]
hkl_planes.extend([pyr_1, pyr_2, pyr_3, pyr_4, pyr_5, pyr_6])
#import sys
#sys.exit(0)

# assembly to gather all parts
assembly = vtk.vtkAssembly()
origin = (a / 2, -a * sqrt(3) / 2., c / 2)

grid = hexagonal_lattice_grid(l, origin=(0., 0., 0.))  # we will use the actor origin
Edges = vtk.vtkExtractEdges()
if vtk.vtkVersion().GetVTKMajorVersion() > 5:
    Edges.SetInputData(grid)
else:
    Edges.SetInput(grid)
Tubes = vtk.vtkTubeFilter()
Tubes.SetInputConnection(Edges.GetOutputPort())
Tubes.SetRadius(0.02 * a)
Tubes.SetNumberOfSides(6)
Tubes.UseDefaultNormalOn()
Tubes.SetDefaultNormal(.577, .577, .577)
# Create the mapper and actor to display the cell edges.
TubeMapper = vtk.vtkPolyDataMapper()
TubeMapper.SetInputConnection(Tubes.GetOutputPort())
hcp_edges = vtk.vtkActor()
hcp_edges.SetMapper(TubeMapper)
assembly.AddPart(hcp_edges)

# add slip planes (with normal)
for hkl in hkl_planes:
    plane = vtk.vtkPlane()
    plane.SetOrigin(origin)
    plane.SetNormal(hkl.normal())
    print('normal is', plane.GetNormal())
    hklplaneActor = add_plane_to_grid(plane, grid, origin, opacity=0.5, show_normal=True)
    assembly.AddPart(hklplaneActor)
s3d.add(assembly)

# apply translation and orientation to the whole assembly
apply_translation_to_actor(assembly, -np.array(origin))


# add axes actor
axes = axes_actor(0.5, fontSize=50)
s3d.add(axes)

# set up camera and render
cam = setup_camera(size=(a, a, c))
center = 1 * (l.matrix[0] + l.matrix[1]) + 0.5 * l.matrix[2]
cam.SetFocalPoint(center)
cam.SetFocalPoint(0., 0., 0.)
cam.SetPosition(4 * a, -2 * a, 2.5 * a)
cam.Dolly(0.9)
s3d.set_camera(cam)
s3d.render()

"""
    PrismaticA
      U   V   T   W | H   K   I   L
      2  -1  -1   0   0   1  -1   0
      1   1  -2   0  -1   1   0   0
      1  -2   1   0  -1   0   1   0

    PyramidalA
      U   V   T   W | H   K   I   L
      2  -1  -1   0   0   1  -1   1
      1   1  -2   0  -1   1   0   1
      1  -2   1   0  -1   0   1  -1
      1  -2   1   0  -1   0   1   1
      1   1  -2   0  -1   1   0  -1
      2  -1  -1   0   0  -1   1   1

    PyramidalCA
       U   V   T   W | H   K   I   L
       2  -1  -1   3  -1   1   0   1
      -1   2  -1  -3  -1   1   0   1
       2  -1  -1  -3  -1   0   1  -1
      -1  -1   2   3  -1   0   1  -1
       2  -1  -1   3  -1   0   1   1
       1   1  -2   3  -1   0   1   1
       2  -1  -1  -3  -1   1   0  -1
       1  -2   1  -3  -1   1   0  -1
       1   1  -2   3   0  -1   1   1
       1  -2   1  -3   0  -1   1   1
      -1  -1   2   3   0   1  -1   1
      -1   2  -1  -3   0   1  -1   1

Also, here are the second order pyramidal.

   U   V   T   W | H   K   I   L
   2  -1  -1   3  -2   1   1   2
   2  -1  -1  -3  -2   1   1  -2
   1   1  -2   3  -1  -1   2   2
   1  -2   1  -3  -1   2  -1  -2
  -1  -1   2   3   1   1  -2   2
  -1   2  -1  -3   1  -2   1  -2

"""