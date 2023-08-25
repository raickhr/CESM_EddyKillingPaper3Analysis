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
            'NLMon45degAndMinus45degStrain_{0:04d}km.nc'.format(
                ell)

        NLM_neg_strain = np.zeros(
            (ndays, ylen, xlen), dtype=float)
        NLM_pos_strain = np.zeros(
            (ndays, ylen, xlen), dtype=float)

        NLMfileName = 'NLmodelEP_{0:04d}km.nc'.format(ell)
        ds_2 = Dataset(ellFold + NLMfileName)

        # slopeFile = 'slopeAndCorr2D_{0:04d}km.nc'.format(ell)
        # slopeDS = Dataset(ellFold + slopeFile)
        # slope2d = np.array(slopeDS.variables['slope'])

        fileName = 'okuboWeissAndStrainDirec_{0:04d}km_AllDay.nc'.format(
            ell)

        ds = Dataset(ellFold + fileName)

        for dayIDX in range(ndays):

            print('working in day {0:d} for ell {1:d}km'.format(
                dayIDX + 1, ell))

            # okuboWeiss_vel = np.array(ds.variables['okuboWeiss_vel'][dayIDX, :, :])

            theta1_vel = np.array(ds.variables['theta1_vel'][dayIDX, :, :])

            # strainDom_vel = okuboWeiss_vel > 0  # -2e-10

            tol = 0
            strainPos = theta1_vel > 0 + tol
            strainNeg = ~strainPos

            # negThetaMask = np.logical_and(strainNeg, strainDom_vel)
            # posThetaMask = np.logical_and(strainPos, strainDom_vel)

            negThetaMask = strainNeg
            posThetaMask = strainPos

            ###############################################################################################################

            NLM_str = np.array(
                ds_2.variables['NLmodel_EPCg_strain'][dayIDX, :, :]).copy()

            NLM_neg_strain[dayIDX, :, :] = NLM_str.copy()
            NLM_pos_strain[dayIDX, :, :] = NLM_str.copy()

            NLM_pos_strain[dayIDX, ~posThetaMask] = float('nan')
            NLM_neg_strain[dayIDX, ~negThetaMask] = float('nan')

        wds = Dataset(writeFileName, 'w', format='NETCDF4')

        wds.createDimension('Y', ylen)
        wds.createDimension('X', xlen)
        wds.createDimension('time', None)

        cdftime = wds.createVariable('time', float, ('time'))
        cdftime.units = 'days since 0050-10-01 00:00:00'
        timeArr = np.arange(ndays)*7+3
        cdftime[:] = timeArr[:]


        createVariableAndWrite(wds, 'negThetaNLMstr', 'ergs/sec/m^2',
                               'NLM_str on tau theta range 0 to -90 deg',
                               ('time', 'Y', 'X'),
                               float,
                               NLM_neg_strain)

        createVariableAndWrite(wds, 'posThetaNLMstr', 'ergs/sec/m^2',
                               'NLM_str on tau theta outside range 0 to 90 deg',
                               ('time', 'Y', 'X'),
                               float,
                               NLM_pos_strain)

        wds.close()

