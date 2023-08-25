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

    maskNames = list(ds_masks.variables.keys())
    maskNames = maskNames[2:len(maskNames)]
    nregions = len(maskNames)

    masks = np.zeros((nregions, ylen, xlen), dtype=float)

    for i in range(nregions):
        masks[i, :, :] = np.array(ds_masks.variables[maskNames[i]])

    globalRegion = masks[maskNames.index('Global'), :, :]
    landMask = KMT < 1

    maskNames = np.array(maskNames, dtype='object')

    for ellIDX in range(nell):
        ell = ellList[ellIDX]
        ellFold = rootFolder + 'test_7dayAvg_NanLand_test/Output/{0:d}km_smth/'.format(ell)

        writeFileName = ellFold + \
            'NLMonPosAndNegVort_{0:04d}km.nc'.format(
                ell)

        NLM_pos_vort = np.zeros((ndays, ylen, xlen), dtype=float)
        NLM_neg_vort = np.zeros((ndays, ylen, xlen), dtype=float)

        NLMfileName = 'NLmodelEP_{0:04d}km.nc'.format(ell)
        ds_2 = Dataset(ellFold + NLMfileName)

        # slopeFile = 'slopeAndCorr2D_{0:04d}km.nc'.format(ell)
        # slopeDS = Dataset(ellFold + slopeFile)
        # slope2d = np.array(slopeDS.variables['slope'])

        fileName = 'okuboWeissAndStrainDirec_{0:04d}km_AllDay.nc'.format(
            ell)

        ds = Dataset(ellFold + fileName)

        for dayIDX in range(ndays):
            # fileName = 'okuboWeissAndStrainDirec_{0:04d}km_day{1:03d}.nc'.format(
            #     ell, dayIDX+1)

            print('working in day {0:d} for ell {1:d}km'.format(
                dayIDX+1, ell))

            # okuboWeiss_vel = np.array(ds.variables['okuboWeiss_vel'][dayIDX, :, :])

            omega_vel = np.array(ds.variables['vorticity'][dayIDX, :, :])

            # vortDom_vel = okuboWeiss_vel <= 0 #2e-10

            # posVortMask = np.logical_and(omega_vel >= 0, vortDom_vel)
            # negVortMask = np.logical_and(omega_vel < 0, vortDom_vel)

            posVortMask = omega_vel > 0
            negVortMask = ~posVortMask

            ###############################################################################################################

            NLM_rot = np.array(
                ds_2.variables['NLmodel_EPCg_rot'][dayIDX, :, :]).copy()

            NLM_pos_vort[dayIDX, :, :] = NLM_rot.copy()
            NLM_neg_vort[dayIDX, :, :] = NLM_rot.copy()

            NLM_pos_vort[dayIDX, ~posVortMask] = float('nan')
            NLM_neg_vort[dayIDX, ~negVortMask] = float('nan')

        wds = Dataset(writeFileName, 'w', format='NETCDF4')

        wds.createDimension('Y', ylen)
        wds.createDimension('X', xlen)
        wds.createDimension('time', None)

        cdftime = wds.createVariable('time', float, ('time'))
        cdftime.units = 'days since 0050-10-01 00:00:00'
        timeArr = np.arange(ndays)*7+3
        cdftime[:] = timeArr[:]


        createVariableAndWrite(wds, 'posVortNLMrot', 'watts/m^2',
                               'NLM_rot on positive Vorticity',
                               ('time', 'Y', 'X'),
                               float,
                               NLM_pos_vort)

        createVariableAndWrite(wds, 'negVortNLMrot', 'watts/m^2',
                               'NLM_rot on negative Vorticity',
                               ('time', 'Y', 'X'),
                               float,
                               NLM_neg_vort)

        # createVariableAndWrite(wds, 'posVortNLMrot_AreaInt', 'watts',
        #                       'NLM_rot on positive Vorticity',
        #                       ('time', 'Y', 'X'),
        #                       float,
        #                       NLM_pos_vort_AreaInt)

        # createVariableAndWrite(wds, 'negVortNLMrot_AreaInt', 'watts',
        #                       'NLM_rot on negative Vorticity',
        #                       ('time', 'Y', ''),
        #                       float,
        #                       NLM_neg_vort_AreaInt)

        wds.close()

