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
    rootFolder = '/pscratch/sd/s/srai/nccs/runGeos7dayAvg/test_7dayAvg_NanLand_test/'
    ds_masks = Dataset(rootFolder + 'input/RegionMasks.nc')
    GridFile = rootFolder + '/input/tripoleGridCreated.nc'


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
        ellFold = rootFolder + 'Output/{0:d}km_smth/'.format(ell)

        writeFileName = ellFold + \
            'avgEPonPosAndNegVort_{0:04d}km.nc'.format(
                ell)

        fileName = 'EPonPosAndNegVort_{0:04d}km.nc'.format(
            ell)

        print('working in {0:s} file for ell {1:d}km'.format(
            fileName, ell))

        ds = Dataset(ellFold + fileName)

        avgEP_pos_vort = np.zeros((ylen, xlen), dtype=float)
        avgEP_neg_vort = np.zeros((ylen, xlen), dtype=float)

        avgEP_pos_vort_AreaInt = np.zeros((ylen, xlen), dtype=float)
        avgEP_neg_vort_AreaInt = np.zeros((ylen, xlen), dtype=float)

        for day in range(ndays):
            print('reading day', day)
            EP_pos_vort = np.array(ds.variables['posVort_EP'][day, :, :])
            EP_neg_vort = np.array(ds.variables['negVort_EP'][day, :, :])

            EP_pos_vort[np.isnan(EP_pos_vort)] = 0.0
            EP_neg_vort[np.isnan(EP_neg_vort)] = 0.0

            avgEP_pos_vort += EP_pos_vort
            avgEP_neg_vort += EP_neg_vort

            # np.nanmean(EP_pos_vort_AreaInt, axis=0)
            avgEP_pos_vort_AreaInt += avgEP_pos_vort * UAREA
            # np.nanmean(EP_neg_vort_AreaInt, axis=0)
            avgEP_neg_vort_AreaInt += avgEP_neg_vort * UAREA

        avgEP_pos_vort /= ndays
        avgEP_neg_vort /= ndays
        avgEP_pos_vort_AreaInt /= ndays
        avgEP_neg_vort_AreaInt /= ndays

        wds = Dataset(writeFileName, 'w', format='NETCDF4')

        wds.createDimension('Y', ylen)
        wds.createDimension('X', xlen)
        # wds.createDimension('time', None)

        createVariableAndWrite(wds, 'posVort_EP', 'ergs/cm^2/sec',
                               'EP on positive Vorticity',
                               ('Y', 'X'),
                               float,
                               avgEP_pos_vort)

        createVariableAndWrite(wds, 'negVort_EP', 'ergs/cm^2/sec',
                               'EP on negative Vorticity',
                               ('Y', 'X'),
                               float,
                               avgEP_neg_vort)

        createVariableAndWrite(wds, 'posVort_EP_AreaInt', 'ergs/sec',
                               'EP on positive Vorticity',
                               ('Y', 'X'),
                               float,
                               avgEP_pos_vort_AreaInt)

        createVariableAndWrite(wds, 'negVort_EP_AreaInt', 'ergs/sec',
                               'EP on negative Vorticity',
                               ('Y', 'X'),
                               float,
                               avgEP_neg_vort_AreaInt)

        wds.close()


