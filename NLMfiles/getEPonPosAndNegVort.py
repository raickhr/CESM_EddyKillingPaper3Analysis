from netCDF4 import Dataset
import numpy as np


def createVariableAndWrite(wds, varName, varUnits, varLongName, varDimension, varType, varArr):
    cdfVar = wds.createVariable(varName, varType, varDimension)
    cdfVar.units = varUnits
    cdfVar.long_name = varLongName
    cdfVar[:] = varArr[:]


def setMaskZeroToNan(ds, varName):
    var = np.array(ds.variables[varName][0, :, :]).copy()
    var[var == 0] = float('nan')
    return var


if __name__ == '__main__':
    #rootFolder = '/discover/nobackup/projects/mesocean/srai/runGeos7dayAvg/'

    #ds_masks = Dataset(rootFolder + 'tavgInput/RegionMasks.nc')
    #GridFile = rootFolder + 'tavgInput/newGlobalGrid_tripolePOP_0.1deg.nc'
    
    rootFolder = '/pscratch/sd/s/srai/nccs/runGeos7dayAvg/'
    ds_masks = Dataset(rootFolder + 'test_7dayAvg_NanLand_test/input/RegionMasks.nc')
    GridFile = rootFolder + 'test_7dayAvg_NanLand_test/input/tripoleGridCreated.nc'

    gridDS = Dataset(GridFile)
    UAREA = np.array(gridDS.variables['UAREA'])
    KMT = np.array(gridDS.variables['KMT'])
    ylen, xlen = np.shape(UAREA)
    DX = np.array(gridDS.variables['DXU'])
    DY = np.array(gridDS.variables['DYU'])

    # ellList = np.arange(60, 450, 20)
    # ellList = np.append(np.array([50]), ellList, axis=0)
    # ellList = np.append(ellList, np.arange(500, 1050, 100), axis=0)
    ellList = np.array([100])
    nell = len(ellList)
    ndays = 365

    landMask = KMT < 1

    for ellIDX in range(nell):
        ell = ellList[ellIDX]
        ellFold = rootFolder + 'test_7dayAvg_NanLand_test/Output/{0:d}km_smth/'.format(ell)

        writeFileName = ellFold + \
            'EPonPosAndNegVort_{0:04d}km.nc'.format(
                ell)

        EP_pos_vort = np.zeros((ndays, ylen, xlen), dtype=float)
        EP_neg_vort = np.zeros((ndays, ylen, xlen), dtype=float)

        fileName = 'filtered_{0:04d}.nc'.format(ell)
        ds_2 = Dataset(ellFold +fileName)

        okuboFileName = 'okuboWeissAndStrainDirec_{0:04d}km_AllDay.nc'.format(
            ell)

        ds = Dataset(ellFold + okuboFileName)

        for dayIDX in range(ndays):
            # fileName = 'okuboWeissAndStrainDirec_{0:04d}km_day{1:03d}.nc'.format(
            #     ell, dayIDX+1)

            print('working in day {0:d} for ell {1:d}km'.format(
                dayIDX+1, ell))

            okuboWeiss_vel = np.array(
                ds.variables['okuboWeiss_vel'][dayIDX, :, :])

            omega_vel = np.array(ds.variables['vorticity'][dayIDX, :, :])

            vortDom_vel = okuboWeiss_vel <= 0  # 2e-10

            posVortMask = np.logical_and(omega_vel >= 0, vortDom_vel)
            negVortMask = np.logical_and(omega_vel < 0, vortDom_vel)

            ###############################################################################################################

            EP_rot = np.array(
                ds_2.variables['EddyPowerPerArea'][dayIDX, :, :]).copy()

            EP_pos_vort[dayIDX, :, :] = EP_rot.copy()
            EP_neg_vort[dayIDX, :, :] = EP_rot.copy()

            EP_pos_vort[dayIDX, ~posVortMask] = float('nan')  # 0.0
            EP_neg_vort[dayIDX, ~negVortMask] = float('nan')  # 0.0

        wds = Dataset(writeFileName, 'w', format='NETCDF4')

        wds.createDimension('Y', ylen)
        wds.createDimension('X', xlen)
        wds.createDimension('time', None)

        cdftime = wds.createVariable('time', float, ('time'))
        cdftime.units = 'days since 0050-01-01 00:00:00'
        timeArr = np.arange(ndays)*7+3
        cdftime[:] = timeArr[:]

        createVariableAndWrite(wds, 'posVort_EP', 'ergs/sec/cm^2',
                               'EP_rot on positive Vorticity',
                               ('time', 'Y', 'X'),
                               float,
                               EP_pos_vort)

        createVariableAndWrite(wds, 'negVort_EP', 'ergs/sec/cm^2',
                               'EP_rot on negative Vorticity',
                               ('time', 'Y', 'X'),
                               float,
                               EP_neg_vort)

        wds.close()

