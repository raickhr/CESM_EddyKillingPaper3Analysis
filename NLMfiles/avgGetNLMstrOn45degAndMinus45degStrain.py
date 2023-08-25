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
            'avgNLMon45degAndMinus45degStrain_{0:04d}km.nc'.format(
                ell)

        fileName = 'NLMon45degAndMinus45degStrain_{0:04d}km.nc'.format(
            ell)

        print('working in {0:s} file for ell {1:d}km'.format(
            fileName, ell))

        ds = Dataset(ellFold + fileName)

        avgNLM_pos_strain = np.zeros((ylen, xlen), dtype=float)
        avgNLM_neg_strain = np.zeros((ylen, xlen), dtype=float)

        avgNLM_pos_strain_AreaInt = np.zeros((ylen, xlen), dtype=float)
        avgNLM_neg_strain_AreaInt = np.zeros((ylen, xlen), dtype=float)

        for day in range(ndays):
            print('reading day', day)
            NLM_pos_strain = np.array(ds.variables['posThetaNLMstr'][day, :, :])
            NLM_neg_strain = np.array(ds.variables['negThetaNLMstr'][day, :, :])

            NLM_pos_strain[np.isnan(NLM_pos_strain)] = 0.0
            NLM_neg_strain[np.isnan(NLM_neg_strain)] = 0.0

            avgNLM_pos_strain += NLM_pos_strain
            avgNLM_neg_strain += NLM_neg_strain

            # np.nanmean(NLM_pos_strain_AreaInt, axis=0)
            avgNLM_pos_strain_AreaInt += avgNLM_pos_strain * UAREA
            # np.nanmean(NLM_neg_strain_AreaInt, axis=0)
            avgNLM_neg_strain_AreaInt += avgNLM_neg_strain * UAREA

        avgNLM_pos_strain /= ndays
        avgNLM_neg_strain /= ndays
        avgNLM_pos_strain_AreaInt /= ndays
        avgNLM_neg_strain_AreaInt /= ndays

        wds = Dataset(writeFileName, 'w', format='NETCDF4')

        wds.createDimension('Y', ylen)
        wds.createDimension('X', xlen)
        # wds.createDimension('time', None)

        createVariableAndWrite(wds, 'posThetaNLMstr', 'ergs/cm^2/sec',
                               'NLM_str on thetaVel > 0',
                               ('Y', 'X'),
                               float,
                               avgNLM_pos_strain)

        createVariableAndWrite(wds, 'negThetaNLMstr', 'ergs/cm^2/sec',
                               'NLM_str on thetaVel < 0',
                               ('Y', 'X'),
                               float,
                               avgNLM_neg_strain)

        createVariableAndWrite(wds, 'posThetaNLMstr_AreaInt', 'ergs/sec',
                               'NLM_str on thetaVel > 0',
                               ('Y', 'X'),
                               float,
                               avgNLM_pos_strain_AreaInt)

        createVariableAndWrite(wds, 'negThetaNLMstr_AreaInt', 'ergs/sec',
                               'NLM_str on thetaVel < 0',
                               ('Y', 'X'),
                               float,
                               avgNLM_neg_strain_AreaInt)

        wds.close()

