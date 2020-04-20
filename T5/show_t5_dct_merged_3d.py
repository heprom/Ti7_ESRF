import os
import h5py
import vtk
from pymicro.view.vtk_utils import numpy_support
from pymicro.view.vtk_utils import *
from pymicro.view.scene3d import *
from pymicro.file.file_utils import HST_read
from scipy import ndimage

show_sample = False

# load our data set
ds = 2  # downsample
#data_dir = 'id11/t5_/t5_dct_bot_/5_reconstruction'
#f = h5py.File(os.path.join(data_dir, 'volume_dilated.mat'))
#grain_ids = f['grains'][::ds, 100:590:ds, 100:630:ds].transpose(2, 1, 0)
grain_ids = HST_read('t5_test_cat.raw')[::ds, ::ds, ::ds]
dims = np.array(grain_ids.shape)
grain_ids[:dims[0] // 2, :dims[1] // 2, :] = 0
#f.close()

# create our vtk actors
all_grains = show_grains(grain_ids)
box = box_3d(size=dims)
axes = axes_actor(100, axisLabels=('X', 'Y', 'Z'), fontSize=50)

# create the 3D scene
s3d = Scene3D(display=False, ren_size=(800, 950), name='t5_dct_merged')
s3d.add(all_grains)
s3d.add(box)
s3d.add(axes)

if show_sample:
    mask = HST_read('t5_test_cat_mask.raw')[::ds, ::ds, ::ds]
    sample = show_array(mask)
    sample.GetProperty().SetOpacity(0.3)
    s3d.add(sample)

# render
cam = setup_camera(size=dims)
cam.SetFocalPoint(0.5 * dims)
cam.SetViewUp(0, 0, 1)
cam.SetPosition(-2 * dims[0], -1.5 * dims[0], 2 * dims[2])
s3d.set_camera(cam)
s3d.render()
