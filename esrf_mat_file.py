import numpy as np
import h5py
import scipy.ndimage
import pymicro
import vtk


class esrf_mat_file:
    def __init__(self,file):
        self.file = file # name of file
        print('file = '+self.file)
        self.header = None
        self.headerSize = 512
        self.dataType = np.float64
        self.dimx = 1
        self.dimy = 1
        self.dimz = 0
        self.read_volume()
        # f = h5py.File(self.file)
        # print(list(f.keys()))
        # self.vol = 0
        # self.read_volume
        # self.get_header()
        # self.parse_header()
        # self.esrf_to_numpy_datatype()
        # self.read_image_data()

    def extract_grain(self,grainId, pad = None):
        print('Extracting grain '+str(grainId)+' to self.currentGrain')
        self.grainId = grainId
        grain = -1
        r_vec = -1
        if pad == None:
            pad = 10
        try:
            size = self.vol.shape
            mask = self.vol == grainId
            vol = (mask) * self.vol
            xyz = np.where(mask == 1)
            x_range = np.array([xyz[0].min(), xyz[0].max()])
            y_range = np.array([xyz[1].min(), xyz[1].max()])
            z_range = np.array([xyz[2].min(), xyz[2].max()])

            grain = vol[x_range[0] - pad:x_range[1] + pad, y_range[0] - pad:y_range[1] + pad,
                       z_range[0] - pad:z_range[1] + pad]
            self.currentGrain = grain
        except Exception as err:
            print('Error extract_grain isolating grain...')
            print("error {}".format(err))

        try:
            print('Getting rodrigues vector for grain '+str(self.grainId))
            nm = self.file
            grain_name = nm[:nm.find('5_reconstruction')] + '4_grains/phase_01/grain_' + str(grainId).zfill(4) + '.mat'
            gf = h5py.File(grain_name)
            r_vec = np.copy(gf['R_vector'])
            gf.close()
            self.currentRodriguesVector = r_vec.flatten()
        except Exception as err:
            print('Error in extract_grain getting rodrigues vector...')
            print("error {}".format(err))

    def extract_multiple_grains(self,grainIDs,pad = None):
        print('Extracting grains to self.multipleGrains')
        self.grainIDs = grainIDs
        grains = -1
        r_vecs = {}
        if pad == None:
            pad = 10
        try:
            size = self.vol.shape
            mask = np.isin(self.vol,grainIDs)
            vol = (mask) * self.vol
            xyz = np.where(mask == 1)
            x_range = np.array([xyz[0].min(), xyz[0].max()])
            y_range = np.array([xyz[1].min(), xyz[1].max()])
            z_range = np.array([xyz[2].min(), xyz[2].max()])

            grains = vol[x_range[0] - pad:x_range[1] + pad, y_range[0] - pad:y_range[1] + pad,
                       z_range[0] - pad:z_range[1] + pad]
            self.multipleGrains = grains
        except Exception as err:
            print('Error extract_grain isolating grain...')
            print("error {}".format(err))

    def get_grain_orientations_dictionary(self):
        ids = np.unique(self.vol)
        ids = ids[np.where(ids > 0)]
        r_vecs = {}
        print('Getting rodrigues vector dictionary for '+str(ids.shape)+' grains')
        for ID in ids:
            try:
                # print('Getting rodrigues vector dictionary for grains')
                nm = self.file
                grain_name = nm[:nm.find('5_reconstruction')] + '4_grains/phase_01/grain_' + str(ID).zfill(4) + '.mat'
                gf = h5py.File(grain_name)
                r_vec = np.copy(gf['R_vector'])
                gf.close()
                RodriguesVector = r_vec.flatten()
                r_vecs[ID] = RodriguesVector
            except Exception as err:
                print('Error in extract_grain getting rodrigues vector...')
                print("error {}".format(err))
        self.r_vecs = r_vecs

    def read_volume(self):
        print('reading volume...')
        f = h5py.File(self.file, 'r')
        print(list(f.keys()))
        vol = np.copy(f['vol'])
        self.vol = vol
        f.close()
        self.dimx = vol.shape[0]
        self.dimy = vol.shape[1]
        self.dimz = vol.shape[2]
        print('shape ='+str(self.vol.shape))


    def get_header(self):
        with open(self.file, 'rb') as f:
            head = f.read(512)
        self.header = head.decode("utf-8")
        self.parse_header()

    def parse_header(self):
        head = self.header
        head = head.replace('{', '')
        head = head.replace('\n', '')
        head = head.replace('{','')
        d = list()
        header_dictionary = {}
        for txt in head.split(';'):
            print(txt)
            if '=' in txt:
                print('found =')
                tmp = txt.split('=')
                d.append(tmp)
                d1 = tmp[0].strip()
                d2 = tmp[1].strip()
                header_dictionary[d1.lower()] = d2
        try:
            size = np.float(header_dictionary['size'])
            size_header = os.path.getsize(self.file) - size
            self.size = size
            self.size_header_in_bytes = size_header
            self.headrSize = size_header
        except:
            print('error getting size in bytes')
        try:
            dimx = np.float(header_dictionary['dim_1'])
            print('dimx ='+str(dimx))
            self.dimx = np.int(dimx)
            dimy = np.float(header_dictionary['dim_2'])
            print('dimy ='+str(dimy))
            self.dimy = np.int(dimy)
            dimz = np.float(header_dictionary['dim_3'])
            print('dimz ='+str(dimz))
            self.dimz = np.int(dimz)
        except:
            print('error getting dimensions')
        try:
            dtype = header_dictionary['datatype']
            self.dataTypeString = dtype
        except:
            print('error getting data type')

    def get_microstructure_from_grains(self,ids = None,name = None):
        if ids == None:
            ids = np.unique(self.vol)
            ids = ids[np.where(ids > 0)]

        if name == None:
            name = 'empty'
        micro = pymicro.crystal.microstructure.Microstructure(name = name)
        microMesh = vtk.vtkMultiBlockDataSet()
        microMesh.SetNumberOfBlocks(len(ids))

        for i, id in enumerate(ids):
            print('Working on grain '+str(id))
            mask = self.vol == id
            grain_vol = (mask) * self.vol
            # grain_vol = (self.vol == id).astype(int)
            try:
                print(self.r_vecs[id])
                r = np.copy(self.r_vecs[id])
                orientation = pymicro.crystal.microstructure.Orientation.from_rodrigues(r)
            except:
                print('Orientation not found for grain '+str(id)+', skipping.')
                continue
            grain = pymicro.crystal.microstructure.Grain(id, orientation)
            grain.position = scipy.ndimage.measurements.center_of_mass(grain_vol)
            grain.volume = scipy.ndimage.measurements.sum(grain_vol)
            grain.add_vtk_mesh(grain_vol, contour = False)
            print('Adding mesh for grain '+str(id)+'.')
            microMesh.SetBlock(i,grain.vtkmesh)
            micro.grains.append(grain)
        micro.SetVtkMesh(microMesh)
        self.micro = micro
        # micro.save()


    def get_grain_centroids(self):
        ids = np.unique(self.vol)
        ids = ids[np.where(ids > 0)]
        centers_of_mass = {}
        print('Getting centroid for '+str(ids.shape)+' grains')
        for ID in ids:
            try:
                print('Getting centroid for grain '+str(ID))
                mask = self.vol == ID
                vol = (mask) * self.vol
                ## This is obviously not the best way to do this. should just calculate without
                ## any reordering.
                vol_xyz = np.copy(vol)
                vol_xyz = vol_xyz.swapaxes(0, 2)
                vol_xyz = np.copy(vol_xyz, 'C_CONTIGUOUS')
                com = scipy.ndimage.measurements.center_of_mass(vol_xyz, vol_xyz)
                centers_of_mass[ID] = com
            except Exception as err:
                print('Error in extract_grain getting rodrigues vector...')
                print("error {}".format(err))
        self.centers_of_mass = centers_of_mass

    def esrf_to_numpy_datatype(self):
        print('getting datatype')