import os

# ellList=np.arange(40,380,20)
# ellList=np.append(np.array([50]),ellList, axis=0)
ellList = [100]#,180]
#223.nc_100._Filtered.nc
for ell in ellList:
    for i in range(10,12): #366):
        fileName = '{1:d}km_30km/{0:d}_30km.nc_{1:d}._Filtered.nc'.format(
            i, ell)

        newFileName = '{1:d}km_30km/{0:03d}.nc_{1:04d}_Filtered.nc'.format(
            i,  ell)

        cmd = 'mv {0:s} {1:s}'.format(fileName, newFileName)
        print(cmd)
        os.system(cmd)

