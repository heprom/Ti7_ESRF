import os
import h5py
import vtk
from pymicro.view.vtk_utils import numpy_support
from pymicro.view.vtk_utils import *
from pymicro.view.scene3d import *
from scipy import ndimage

sample = 't5'
data_dir = '/home/proudhon/ucsb/Ti7_ESRF/id11/t5_/t5_dct_cen_1_/5_reconstruction'
f = h5py.File(os.path.join(data_dir, 'phase_01_vol.mat'))
grain_ids = f['vol'].value[::3, 100:540:3, 80:560:3].transpose(2, 1, 0)
dims = np.array(grain_ids.shape)

grains = show_grains(grain_ids)
box = box_3d(size=dims)
axes = axes_actor(50, axisLabels=('X', 'Y', 'Z'), fontSize=50)
s3d = Scene3D(display=True, ren_size=(800, 800), name='%s_dct_cen_3d' % sample)
s3d.add(grains)
#s3d.add(box)
s3d.add(axes)
cam = setup_camera(size=dims)
cam.SetFocalPoint(0.5 * dims)
#cam.SetPosition(-dims[0], -dims[0], 1.1*dims[2])
cam.SetViewUp(0, 1, 0)
cam.SetPosition(-2 * dims[0], 0.5 * dims[0], 0.5 * dims[2])
s3d.set_camera(cam)
s3d.render(key_pressed_callback=True)
