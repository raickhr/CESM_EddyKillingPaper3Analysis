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
            'avgEPon45degAndMinus45Strain_{0:04d}km.nc'.format(
                ell)

        fileName = 'EPon45degAndMinus45Strain_{0:04d}km.nc'.format(
            ell)

        print('working in {0:s} file for ell {1:d}km'.format(
            fileName, ell))

        ds = Dataset(ellFold + fileName)

        avgEP_pos_strain = np.zeros((ylen, xlen), dtype=float)
        avgEP_neg_strain = np.zeros((ylen, xlen), dtype=float)

        avgEP_pos_strain_AreaInt = np.zeros((ylen, xlen), dtype=float)
        avgEP_neg_strain_AreaInt = np.zeros((ylen, xlen), dtype=float)

        for day in range(ndays):
            print('reading day', day)
            EP_pos_strain = np.array(
                ds.variables['posThetaEP_str'][day, :, :])
            EP_neg_strain = np.array(
                ds.variables['negThetaEP_str'][day, :, :])

            EP_pos_strain[np.isnan(EP_pos_strain)] = 0.0
            EP_neg_strain[np.isnan(EP_neg_strain)] = 0.0

            avgEP_pos_strain += EP_pos_strain
            avgEP_neg_strain += EP_neg_strain

            # np.nanmean(EP_pos_strain_AreaInt, axis=0)
            avgEP_pos_strain_AreaInt += avgEP_pos_strain * UAREA
            # np.nanmean(EP_neg_strain_AreaInt, axis=0)
            avgEP_neg_strain_AreaInt += avgEP_neg_strain * UAREA

        avgEP_pos_strain /= ndays
        avgEP_neg_strain /= ndays
        avgEP_pos_strain_AreaInt /= ndays
        avgEP_neg_strain_AreaInt /= ndays

        wds = Dataset(writeFileName, 'w', format='NETCDF4')

        wds.createDimension('Y', ylen)
        wds.createDimension('X', xlen)
        # wds.createDimension('time', None)

        createVariableAndWrite(wds, 'posThetaEP_str', 'ergs/cm^2/sec',
                               'EP_str on thetaVel > 0',
                               ('Y', 'X'),
                               float,
                               avgEP_pos_strain)

        createVariableAndWrite(wds, 'negThetaEP_str', 'ergs/cm^2/sec',
                               'EP_str on thetaVel < 0',
                               ('Y', 'X'),
                               float,
                               avgEP_neg_strain)

        createVariableAndWrite(wds, 'posThetaEP_str_AreaInt', 'ergs/sec',
                               'EP_str on thetaVel > 0',
                               ('Y', 'X'),
                               float,
                               avgEP_pos_strain_AreaInt)

        createVariableAndWrite(wds, 'negThetaEP_str_AreaInt', 'ergs/sec',
                               'EP_str on thetaVel < 0',
                               ('Y', 'X'),
                               float,
                               avgEP_neg_strain_AreaInt)

        wds.close()


