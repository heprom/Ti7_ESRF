from math import sqrt
import numpy as np
from pymicro.view.scene3d import Scene3D
from pymicro.view.vtk_utils import *
from vtk.util.colors import white, peacock, tomato, red, green, yellow, black, cyan, magenta
from pymicro.crystal.lattice import HklPlane
from pymicro.crystal.microstructure import Orientation

'''
Create a 3d scene with a series of hexagonal lattices representing different grains.
Hkl planes are added to the lattice and cut by the Z=0 plane to display slip traces.

'''
# Create the 3D scene
base_name = os.path.splitext(__file__)[0]
s3d = Scene3D(display=True, ren_size=(800, 800), name=base_name)

# hexagonal lattice for Ti7Al
a = 0.2931  # nm
c = 0.4694  # nm
l = Lattice.hexagonal(a, c)

# grain 26
o26 = Orientation.from_euler((296.091, 112.686,  24.046))
#o26 = Orientation.cube()
print(o26)
basal = HklPlane(0, 0, 1, lattice=l)
prism_1 = HklPlane(0, 1, 0, lattice=l)
prism_2 = HklPlane(-1, 1, 0, lattice=l)
prism_3 = HklPlane(-1, 0, 0, lattice=l)
hkl_planes = [basal]
#hkl_planes = [prism_1, prism_2, prism_3]
colors = [black, red, green, yellow, black, cyan, magenta]  # correspond to pyplot 'rgykcmbw'

# assembly to gather all parts
assembly = vtk.vtkAssembly()
origin = (a / 2, -a * sqrt(3) / 2., c / 2)

# intersection plane to show the slip traces
inter_plane = vtk.vtkPlane()
inter_plane.SetOrigin(origin)
n_inter = np.array([1, 0, 0])  # to plot slip traces on XY
inter_plane.SetNormal(o26.orientation_matrix().dot(n_inter))  # here the plane is rotated to crystal CS

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
for i, hkl in enumerate(hkl_planes):
    plane = vtk.vtkPlane()
    plane.SetOrigin(origin)
    plane.SetNormal(hkl.normal())
    print('normal is', plane.GetNormal())
    hklplaneActor = add_plane_to_grid(plane, grid, origin, opacity=0.5, show_normal=True)
    # get a reference to the vtkPolyData representing the hkl plane
    hklplanePolyData = hklplaneActor.GetParts().GetItemAsObject(0).GetMapper().GetInput()
    assembly.AddPart(hklplaneActor)

    # cut the rotated plane with the vertical plane to display the trace
    slipTrace = vtk.vtkCutter()
    slipTrace.SetInputData(hklplanePolyData)
    slipTrace.SetCutFunction(inter_plane)
    slipTrace.Update()  # this is a vtkPolyData
    slipTraceMapper = vtk.vtkPolyDataMapper()
    slipTraceMapper.SetInputConnection(slipTrace.GetOutputPort())
    slipTraceActor = vtk.vtkActor()
    slipTraceActor.SetMapper(slipTraceMapper)
    slipTraceActor.GetProperty().SetColor(colors[i])
    slipTraceActor.GetProperty().SetLineWidth(5)
    assembly.AddPart(slipTraceActor)

apply_orientation_to_actor(assembly, o26)
s3d.add(assembly)

# add axes actor
axes = axes_actor(0.5, fontSize=50)
s3d.add(axes)

# set up camera and render
cam = setup_camera(size=(a, a, c))
center = 1 * (l.matrix[0] + l.matrix[1]) + 0.5 * l.matrix[2]
cam.SetFocalPoint(center)
cam.SetFocalPoint(0., 0., 0.)
cam.SetPosition(-5 * a, 0, 0)
cam.SetViewUp(0, -1, 0)
cam.Dolly(0.5)
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