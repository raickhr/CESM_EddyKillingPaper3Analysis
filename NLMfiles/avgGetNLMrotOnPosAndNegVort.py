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
            'avgNLMonPosAndNegVort_{0:04d}km.nc'.format(
                ell)

        fileName = 'NLMonPosAndNegVort_{0:04d}km.nc'.format(
            ell)

        print('working in {0:s} file for ell {1:d}km'.format(
            fileName, ell))

        ds = Dataset(ellFold + fileName)

        avgNLM_pos_vort = np.zeros((ylen, xlen), dtype=float)
        avgNLM_neg_vort = np.zeros((ylen, xlen), dtype=float)

        avgNLM_pos_vort_AreaInt = np.zeros((ylen, xlen), dtype=float)
        avgNLM_neg_vort_AreaInt = np.zeros((ylen, xlen), dtype=float)

        for day in range(ndays):
            print('reading day', day)
            NLM_pos_vort = np.array(ds.variables['posVortNLMrot'][day, :, :])
            NLM_neg_vort = np.array(ds.variables['negVortNLMrot'][day, :, :])

            NLM_pos_vort[np.isnan(NLM_pos_vort)] = 0.0
            NLM_neg_vort[np.isnan(NLM_neg_vort)] = 0.0

            avgNLM_pos_vort += NLM_pos_vort
            avgNLM_neg_vort += NLM_neg_vort

            # np.nanmean(NLM_pos_vort_AreaInt, axis=0)
            avgNLM_pos_vort_AreaInt += avgNLM_pos_vort * UAREA
            # np.nanmean(NLM_neg_vort_AreaInt, axis=0)
            avgNLM_neg_vort_AreaInt += avgNLM_neg_vort * UAREA

        avgNLM_pos_vort /= ndays
        avgNLM_neg_vort /= ndays
        avgNLM_pos_vort_AreaInt /= ndays
        avgNLM_neg_vort_AreaInt /= ndays

        wds = Dataset(writeFileName, 'w', format='NETCDF4')

        wds.createDimension('Y', ylen)
        wds.createDimension('X', xlen)
        # wds.createDimension('time', None)

        createVariableAndWrite(wds, 'posVortNLMrot', 'ergs/cm^2/sec',
                               'NLM_rot on positive Vorticity',
                               ('Y', 'X'),
                               float,
                               avgNLM_pos_vort)

        createVariableAndWrite(wds, 'negVortNLMrot', 'ergs/cm^2/sec',
                               'NLM_rot on negative Vorticity',
                               ('Y', 'X'),
                               float,
                               avgNLM_neg_vort)

        createVariableAndWrite(wds, 'posVortNLMrot_AreaInt', 'ergs/sec',
                               'NLM_rot on positive Vorticity',
                               ('Y', 'X'),
                               float,
                               avgNLM_pos_vort_AreaInt)

        createVariableAndWrite(wds, 'negVortNLMrot_AreaInt', 'ergs/sec',
                               'NLM_rot on negative Vorticity',
                               ('Y', 'X'),
                               float,
                               avgNLM_neg_vort_AreaInt)

        wds.close()

