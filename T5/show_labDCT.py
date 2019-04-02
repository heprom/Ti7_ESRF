import os
import h5py
import vtk
from pymicro.view.vtk_utils import numpy_support
from pymicro.view.vtk_utils import *
from pymicro.view.scene3d import *
from scipy import ndimage

sample = 'T5'
data_dir = '/home/proudhon/data/tomo/rawdata/2018_ma3921/t5_'
f = h5py.File(os.path.join(data_dir, 'Ti7Al_T5_stitched.h5'))
FeatureIds = f['LabDCT/Data/GrainId'].value.transpose(1, 2, 0)
FeatureIds = ndimage.rotate(FeatureIds, 11, axes=(1, 0), order=0, reshape=False)[70:185, 50:165, :]
'''
from matplotlib import pyplot as plt
plt.imshow(FeatureIds[:, :, 100].T, interpolation='nearest')
plt.show()
'''
dims = FeatureIds.shape
IPFColors = (f['LabDCT/Data/IPF001'].value * 255).astype(np.uint8).transpose(1, 2, 0, 3)
IPFColors = ndimage.rotate(IPFColors, 11, axes=(1, 0), order=0, reshape=False)[70:185, 50:165, :, :]
print(IPFColors[100, 100, 100, :])
print('IPFColors.dtype = %s' % IPFColors.dtype)
print(numpy_support.get_vtk_array_type(IPFColors.dtype))
print(IPFColors.shape)

# get lattice parameters
a, b, c, alpha, beta, gamma = f['PhaseInfo/Phase01/UnitCell'].value
from pymicro.crystal.lattice import Lattice
Ti7 = Lattice.hexagonal(a, c)

ipf_grid = numpy_array_to_vtk_grid(IPFColors, cell_data=True)
mask = FeatureIds > 0
sample_mask_array = numpy_support.numpy_to_vtk(np.ravel(mask, order='F').astype(np.uint8), deep=1)
ipf_grid.SetCellVisibilityArray(sample_mask_array)
bounds = ipf_grid.GetBounds()
extract = vtk.vtkExtractGeometry()
extract.SetInputData(ipf_grid)
extract.ExtractInsideOn()
extract.ExtractBoundaryCellsOn()
bbox = vtk.vtkBox()
bbox.SetXMin(bounds[0::2])
bbox.SetXMax(bounds[1::2])
extract.SetImplicitFunction(bbox)
extract.Update()

mapper = vtk.vtkDataSetMapper()
mapper.ScalarVisibilityOn()
mapper.SetScalarModeToUseCellData()
mapper.SetInputConnection(extract.GetOutputPort())
#mapper.SetInputData(ipf_grid)
mapper.Update()
grains = vtk.vtkActor()
grains.SetMapper(mapper)

#gids_grid = numpy_array_to_vtk_grid(FeatureIds, cell_data=True, array_name='FeatureIds')
#gb = show_boundaries(gids_grid, array_name='FeatureIds', write=True)
gb_reader = vtk.vtkXMLPolyDataReader()
gb_reader.SetFileName('gb.vtp')
gb_reader.Update()
gb = edges_actor(gb_reader.GetOutput(), linewidth=4.0, linecolor=black)
  
box = box_3d(size=IPFColors.shape)
axes = axes_actor(30, axisLabels=('X', 'Y', 'Z'), fontSize=50)
s3d = Scene3D(display=False, ren_size=(500, 700), name='%s_3d' % sample)
s3d.add(grains)
s3d.add(gb)
#s3d.add(box)
#s3d.add(axes)
cam = setup_camera(size=dims)
cam.SetFocalPoint(0.5 * dims[0], 0.5 * dims[1], 0.47 * dims[2])
cam.SetPosition(-dims[0], -dims[0], 1.1*dims[2])
cam.Dolly(0.55)
s3d.set_camera(cam)
s3d.render(key_pressed_callback=False)
