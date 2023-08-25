from netCDF4 import Dataset
import numpy as np


ellList = np.arange(20,450,20)
ellList2 = np.arange(500,1050,100)
ellList = np.concatenate((ellList, ellList2))
del ellList2

ylen, xlen = 2400, 3600
nell = len(ellList)

fldLoc = '/discover/nobackup/projects/mesocean/srai/runGeos7dayAvg/test_7dayAvg_NanLand_test/Output/'
wfileName = 'KEallEll_day10.nc'

KE = np.zeros((nell, ylen, xlen), dtype=float)

for ellIdx  in range(nell):
    ell = ellList[ellIdx]
    fileName = fldLoc + '10_smth.nc_{0:d}._Filtered.nc'.format(ell)
    ds = Dataset(fileName)
    u = np.array(ds.variables['UVEL'][:,:])
    v = np.array(ds.variables['VVEL'][:,:])
    KE[ellIdx,:,:] = 0.5*(u**2 +v**2)
    ds.close()

ds = Dataset(wfileName,'w', format='NETCDF4')

ds.createDimension('ell', nell)
ds.createDimension('Y', ylen)
ds.createDimension('X', xlen)

wdsEll = ds.createVariable('ell', float, ('ell',))
wdsEll[:] = ellList[:]

wdsKE = ds.createVariable('KE',float, ('ell','Y','X'))
wdsKE[:,:,:] = KE[:,:,:]

ds.close()
