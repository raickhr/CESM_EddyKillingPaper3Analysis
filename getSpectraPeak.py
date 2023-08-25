from netCDF4 import Dataset 
import numpy as np


ellList = np.arange(20,350,20)
nell = len(ellList)
ylen, xlen = 2400,3600
data = np.zeros((nell, ylen, xlen), dtype=float)
spectra = np.zeros((nell-1, ylen, xlen), dtype=float)
rootDir = '/discover/nobackup/projects/mesocean/srai/runGeos7dayAvg/test_7dayAvg_NanLand_test/Output/'
wfileName = rootDir + 'spectraPeak.nc'

for ellIndx in range(nell):
    ell = ellList[ellIndx]
    ellDir = '{0:d}km/'.format(ell)
    fileName = rootDir + ellDir + '010.nc_{0:04d}_Filtered.nc'.format(ell)
    ds = Dataset(fileName)
    uvel = np.array(ds.variables['UVEL'])
    vvel = np.array(ds.variables['VVEL'])
    if len(uvel.shape) == 3:
        uvel = uvel[0,:,:]
        vvel = vvel[0,:,:]

    data[ellIndx,:,:] = 0.5*(uvel**2 + vvel**2)
    ds.close()

ellArray = np.array(ellList)
midell = (np.roll(ellArray, -1)[0:-1] + ellArray[0:-1])/2
dell = np.roll(ellArray, -1)[0:-1] - ellArray[0:-1]
dKE = np.roll(data, -1)[0:-1] - data[0:-1]
for i in range(nell-1):
    spectra[i,:,:] = -1*midell[i]**2 * dKE[i,:,:]/dell[i]
indicesList = np.argmax(spectra, axis = 0)
spectraPeak = midell[indicesList]

wds = Dataset(wfileName, 'w', format= 'NETCDF4')

wds.createDimension('Y', ylen)
wds.createDimension('X', xlen)

var = wds.createVariable('spectraPeak', float, ('Y','X'))
var[:,:] = spectraPeak[:,:]
wds.close()


    


