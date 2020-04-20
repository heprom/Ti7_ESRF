import os
import h5py
import vtk
from pymicro.view.vtk_utils import numpy_support
from pymicro.view.vtk_utils import *
from pymicro.view.scene3d import *
from scipy import ndimage


# load our data set
ds = 4  # downsample
data_dir = 'id11/t5_/t5_dct_bot_/5_reconstruction'
f = h5py.File(os.path.join(data_dir, 'volume_dilated.mat'))
grain_ids = f['grains'][::ds, 100:590:ds, 100:630:ds].transpose(2, 1, 0)
dims = np.array(grain_ids.shape)
f.close()

# create our vtk actors
sample = show_array(grain_ids)
sample.GetProperty().SetOpacity(0.3)
selected_grains_data = np.zeros_like(grain_ids)
gids = [5, 19]
for gid in gids:
    selected_grains_data += gid * (grain_ids == gid)
selected_grains = show_grains(selected_grains_data, num_colors=167)
all_grains = show_grains(grain_ids)
box = box_3d(size=dims)
axes = axes_actor(50, axisLabels=('X', 'Y', 'Z'), fontSize=50)

# create the 3D scene
s3d = Scene3D(display=False, ren_size=(800, 800), name='t5_dct_bot_3d_1')
#s3d.add(sample)
#s3d.add(selected_grains)
s3d.add(all_grains)
s3d.add(box)
s3d.add(axes)

# render
cam = setup_camera(size=dims)
cam.SetFocalPoint(0.5 * dims)
cam.SetViewUp(0, 0, 1)
cam.SetPosition(-2 * dims[0], -1.5 * dims[0], 2 * dims[2])
s3d.set_camera(cam)
s3d.render()
