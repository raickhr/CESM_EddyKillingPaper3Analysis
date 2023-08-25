import os
import numpy as np
# ellList=np.arange(40,380,20)
# ellList=np.append(np.array([50]),ellList, axis=0)
#ellList = np.arange(20,350, 20)#[20, 30, 50, 80, 120, 150, 180]
#223.nc_100._Filtered.nc
#ellList = np.arange(20,450,20)
#ellList = np.concatenate((ellList, np.arange(500,1030,100)))
ellList = [100]
for ell in ellList:
    for i in range(1, 366):
        fileName = '{1:d}km_smth/{0:d}_smth.nc_{1:d}._Filtered.nc'.format(
            i, ell)

        newFileName = '{1:d}km_smth/{0:03d}.nc_{1:04d}_Filtered.nc'.format(
            i,  ell)

        cmd = 'mv {0:s} {1:s}'.format(fileName, newFileName)
        print(cmd)
        os.system(cmd)

