from netCDF4 import Dataset
import numpy as np


def createVariableAndWrite(wds, varName, varUnits, varLongName, varDimension, varType, varArr):
    cdfVar = wds.createVariable(varName, varType, varDimension)
    cdfVar.units = varUnits
    cdfVar.long_name = varLongName
    cdfVar[:] = varArr[:]


def getAreaAvgMaskedVal(posThetaEP_str,
                        negThetaEP_str,
                        posVort_EP,
                        negVort_EP,
                        UAREA, mask):

    maskedArea = np.nansum(UAREA[mask])

    areaAvg_posThetaEP_str = np.nansum(
        posThetaEP_str[mask] * UAREA[mask])/maskedArea

    areaAvg_negThetaEP_str = np.nansum(
        negThetaEP_str[mask] * UAREA[mask])/maskedArea

    areaAvg_posVort_EP = np.nansum(
        posVort_EP[mask] * UAREA[mask])/maskedArea

    areaAvg_negVort_EP = np.nansum(
        negVort_EP[mask] * UAREA[mask])/maskedArea

    returnArr = np.array([areaAvg_posThetaEP_str,
                          areaAvg_negThetaEP_str,
                          areaAvg_posVort_EP,
                          areaAvg_negVort_EP], dtype=float)

    return returnArr


if __name__ == '__main__':
    rootFolder = '/discover/nobackup/projects/mesocean/srai/runGeos7dayAvg/'

    ds_masks = Dataset(rootFolder + 'tavgInput/RegionMasks.nc')
    GridFile = rootFolder + 'tavgInput/newGlobalGrid_tripolePOP_0.1deg.nc'

    gridDS = Dataset(GridFile)
    UAREA = np.array(gridDS.variables['UAREA'])
    ULAT = np.array(gridDS.variables['ULAT'])
    ULONG = np.array(gridDS.variables['ULONG'])
    KMT = np.array(gridDS.variables['KMT'])
    ylen, xlen = np.shape(UAREA)
    DX = np.array(gridDS.variables['DXU'])
    DY = np.array(gridDS.variables['DYU'])

    mask_20_40N = np.logical_and(ULAT > 20, ULAT < 40)
    mask_5S_5N = np.logical_and(ULAT > -5, ULAT < 5)
    mask_15_40S = np.logical_and(ULAT > -40, ULAT < -15)
    landMask = KMT < 30
    highKE_mask = np.array(ds_masks['HighEKE'][:, :])

    landMaskWithHiKE = np.logical_or(landMask, highKE_mask)
    UAREA[landMaskWithHiKE] = float('nan')

    # ellList = np.arange(60, 450, 20)
    # ellList = np.append(np.array([50]), ellList, axis=0)
    # ellList = np.append(ellList, np.arange(500, 1050, 100), axis=0)
    ellList = np.array([100])
    nell = len(ellList)
    ndays = 365

    for ellIDX in range(nell):
        ell = ellList[ellIDX]
        ellFold = rootFolder + 'tavgOutput/{0:d}km/'.format(ell)

        writeFileName = ellFold + \
            'EP_rotAndstr_LatBandSeries_{0:04d}km.nc'.format(
                ell)

        fileName1 = 'EPonPosAndNegVort_{0:04d}km.nc'.format(
            ell)
        fileName2 = 'EP_on45degAndMinus45Strain_{0:04d}km.nc'.format(
            ell)

        print('working in file {0:s} and {1:s} for ell {2:d}km'.format(
            fileName1, fileName2, ell))

        ds1 = Dataset(ellFold + fileName1)
        ds2 = Dataset(ellFold + fileName2)

        NLMstrNrot_20_40N = np.zeros((ndays, 4), dtype=float)
        NLMstrNrot_15_40S = np.zeros((ndays, 4), dtype=float)
        NLMstrNrot_5S_5N = np.zeros((ndays, 4), dtype=float)

        for day in range(ndays):
            print('reading day', day)
            posVort_EP = np.array(ds1.variables['posVort_EP'][day, :, :])
            negVort_EP = np.array(ds1.variables['negVort_EP'][day, :, :])
            posThetaEP_str = np.array(
                ds2.variables['negThetaEP_str'][day, :, :])
            negThetaEP_str = np.array(
                ds2.variables['posThetaEP_str'][day, :, :])

            posVort_EP[landMaskWithHiKE] = float('nan')
            negVort_EP[landMaskWithHiKE] = float('nan')
            posThetaEP_str[landMaskWithHiKE] = float('nan')
            negThetaEP_str[landMaskWithHiKE] = float('nan')

            mask = mask_20_40N.copy()
            NLMstrNrot_20_40N[day, :] = getAreaAvgMaskedVal(posThetaEP_str,
                                                            negThetaEP_str,
                                                            posVort_EP,
                                                            negVort_EP,
                                                            UAREA, mask)

            mask = mask_15_40S.copy()
            NLMstrNrot_15_40S[day, :] = getAreaAvgMaskedVal(posThetaEP_str,
                                                            negThetaEP_str,
                                                            posVort_EP,
                                                            negVort_EP,
                                                            UAREA, mask)

            mask = mask_5S_5N.copy()
            NLMstrNrot_5S_5N[day, :] = getAreaAvgMaskedVal(posThetaEP_str,
                                                           negThetaEP_str,
                                                           posVort_EP,
                                                           negVort_EP,
                                                           UAREA, mask)

        wds = Dataset(writeFileName, 'w', format='NETCDF4')

        wds.createDimension('time', ndays)

        createVariableAndWrite(wds, 'posThetaEP_str_20_40N', 'ergs/cm^2/sec',
                               'NLM_str on thetaVel > 0',
                               ('time'),
                               float,
                               NLMstrNrot_20_40N[:, 0])

        createVariableAndWrite(wds, 'negThetaEP_str_20_40N', 'ergs/cm^2/sec',
                               'NLM_str on thetaVel < 0',
                               ('time'),
                               float,
                               NLMstrNrot_20_40N[:, 1])

        createVariableAndWrite(wds, 'posVort_EP_20_40N', 'ergs/cm^2/sec',
                               'NLM_rot on vorticity > 0',
                               ('time'),
                               float,
                               NLMstrNrot_20_40N[:, 2])

        createVariableAndWrite(wds, 'negVort_EP_20_40N', 'ergs/cm^2/sec',
                               'NLM_rot on vorticity < 0',
                               ('time'),
                               float,
                               NLMstrNrot_20_40N[:, 3])

        createVariableAndWrite(wds, 'posThetaEP_str_5S_5N', 'ergs/cm^2/sec',
                               'NLM_str on thetaVel > 0',
                               ('time'),
                               float,
                               NLMstrNrot_5S_5N[:, 0])

        createVariableAndWrite(wds, 'negThetaEP_str_5S_5N', 'ergs/cm^2/sec',
                               'NLM_str on thetaVel < 0',
                               ('time'),
                               float,
                               NLMstrNrot_5S_5N[:, 1])

        createVariableAndWrite(wds, 'posVort_EP_5S_5N', 'ergs/cm^2/sec',
                               'NLM_rot on vorticity > 0',
                               ('time'),
                               float,
                               NLMstrNrot_5S_5N[:, 2])

        createVariableAndWrite(wds, 'negVort_EP_5S_5N', 'ergs/cm^2/sec',
                               'NLM_rot on vorticity < 0',
                               ('time'),
                               float,
                               NLMstrNrot_5S_5N[:, 3])

        createVariableAndWrite(wds, 'posThetaEP_str_15_40S', 'ergs/cm^2/sec',
                               'NLM_str on thetaVel > 0',
                               ('time'),
                               float,
                               NLMstrNrot_15_40S[:, 0])

        createVariableAndWrite(wds, 'negThetaEP_str_15_40S', 'ergs/cm^2/sec',
                               'NLM_str on thetaVel < 0',
                               ('time'),
                               float,
                               NLMstrNrot_15_40S[:, 1])

        createVariableAndWrite(wds, 'posVort_EP_15_40S', 'ergs/cm^2/sec',
                               'NLM_rot on vorticity > 0',
                               ('time'),
                               float,
                               NLMstrNrot_15_40S[:, 2])

        createVariableAndWrite(wds, 'negVort_EP_15_40S', 'ergs/cm^2/sec',
                               'NLM_rot on vorticity < 0',
                               ('time'),
                               float,
                               NLMstrNrot_15_40S[:, 3])

        wds.close()

