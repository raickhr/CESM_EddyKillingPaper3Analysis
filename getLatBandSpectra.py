from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

#ellList = np.arange(20, 350, 20)
ellList = np.arange(20,450,20)
ellList = np.concatenate((ellList, np.arange(500,1030,100)))

nell = len(ellList)

ellArray = np.array(ellList)
midell = (np.roll(ellArray, -1)[0:-1] + ellArray[0:-1])/2 * 1000
dell = (np.roll(ellArray, -1)[0:-1] - ellArray[0:-1])*1000 

ndays = 1
ylen, xlen = 2400, 3600

latbands1 = [(-10,10),(10,20),(20,30),(30,40),(40,50),(50,60),(60,70),(70,80),(80,90)]
latbands2 = [(-10,-20),(-20,-30),(-30,-40),(-40,-50),(-50,-60),(-60,-70),(-70,-80),(-80,-90)]
latbands3 = [(sub[1], sub[0]) for sub in latbands2]
latbands = latbands1 + latbands3
latbands.sort(key=lambda x: x[0])

del latbands1, latbands2, latbands3

nbands = len(latbands)

data = np.zeros((nbands, nell), dtype=float)
spectra = np.zeros((ndays, nbands, nell-1), dtype=float)
#spectraPeakEll = np.zeros((ndays, nbands), dtype=float)

rootDir = '/discover/nobackup/projects/mesocean/srai/runGeos7dayAvg/test_7dayAvg_NanLand_test/Output/'
wfileName = rootDir + 'spectraLatBands.nc'

gridDS= Dataset('/discover/nobackup/projects/mesocean/srai/runGeos7dayAvg/test_7dayAvg_NanLand_test/input/tripoleGridCreated.nc')
UAREA = np.array(gridDS.variables['UAREA']) *1/(10000)
ULAT = np.array(gridDS.variables['ULAT'])
KMU = np.array(gridDS.variables['KMU'])
landMask = KMU < 1

maskList = np.zeros((nbands, ylen, xlen), dtype=bool)
for i in range(nbands):
    maskList[i, :, :] = ~np.logical_and(
        ULAT > latbands[i][0], ULAT < latbands[i][1])
    
    


wds = Dataset(wfileName, 'w', format='NETCDF4')

wds.createDimension('latbands', nbands)
wds.createDimension('time', ndays)
wds.createDimension('ell',nell-1)

cdfLatbandVar = wds.createVariable('latbands', str, ('latbands'))
for i in range(nbands):
    cdfLatbandVar[i] = str(latbands[i])

cdfEll = wds.createVariable('ell', float, ('ell'))
cdfEll[:] = midell[:]

var = wds.createVariable('spectra', float, ('time', 'latbands', 'ell'))

for day in range(0,1):#range(10,11): #range(1, ndays):
    #print('day',day)
    for ellIndx in range(nell):
        ell = ellList[ellIndx]
        print('ell', ell)
        ellDir = '{0:d}km_smth/'.format(ell)
        fileName = rootDir + ellDir + \
            '{0:03d}.nc_{1:04d}_Filtered.nc'.format(10, ell)
        print('reading filename', fileName)
        ds = Dataset(fileName)
        uvel = np.array(ds.variables['UVEL'])
        vvel = np.array(ds.variables['VVEL'])
        if len(uvel.shape) == 3:
            uvel = uvel[0, :, :]/100
            vvel = vvel[0, :, :]/100
        KE = 0.5*(uvel**2 + vvel**2)
        KEmask = KE == 0
        for latband_indx in range(nbands):
            mask = np.logical_or(maskList[latband_indx, :, :], KEmask)
            mask = np.logical_or(mask, landMask)
            # plt.pcolormesh(mask)
            # plt.colorbar()
            # plt.show()
            sumKE = np.nansum(KE[~mask] * UAREA[~mask])
            sumArea = np.nansum(UAREA[~mask])
            data[latband_indx, ellIndx] = sumKE/sumArea
        ds.close()

    #data -> nbands, nell
    dKE = np.roll(data, -1, axis=1)[:, 0:-1] - data[:, 0:-1]
    for i in range(nell-1):
        # spectra -> ndays, nbands, nell-1
        spectra[day, :, i] = -1*midell[i]**2 * dKE[:, i]/dell[i]
    # indicesList = np.argmax(spectra, axis=0)
    # spectraPeakEll[day, :, :] = midell[indicesList]

    var[:, :, :] = spectra[day, :, :]

wds.close()

