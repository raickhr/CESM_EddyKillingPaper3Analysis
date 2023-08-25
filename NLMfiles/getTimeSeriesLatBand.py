from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt


def createVariableAndWrite(wds, varName, varUnits, varLongName, varDimension, varType, varArr):
    cdfVar = wds.createVariable(varName, varType, varDimension)
    cdfVar.units = varUnits
    cdfVar.long_name = varLongName
    cdfVar[:] = varArr[:]


def getAreaAvgMaskedVal(posThetaNLMstr,
                        negThetaNLMstr,
                        posVortNLMrot,
                        negVortNLMrot,
                        UAREA, mask):
    
    maskedArea = np.nansum(UAREA[mask])
    
    areaAvg_posThetaNLMstr = np.nansum(
        posThetaNLMstr[mask] * UAREA[mask])/maskedArea
    
    areaAvg_negThetaNLMstr = np.nansum(
        negThetaNLMstr[mask] * UAREA[mask])/maskedArea
    
    areaAvg_posVortNLMrot = np.nansum(
        posVortNLMrot[mask] * UAREA[mask])/maskedArea
    
    areaAvg_negVortNLMrot = np.nansum(
        negVortNLMrot[mask] * UAREA[mask])/maskedArea
    
    returnArr = np.array([areaAvg_posThetaNLMstr,
                          areaAvg_negThetaNLMstr,
                          areaAvg_posVortNLMrot,
                          areaAvg_negVortNLMrot], dtype=float)

    #returnArr = np.stack([areaAvg_posThetaNLMstr, 
    #                      areaAvg_negThetaNLMstr, 
    #                      areaAvg_posVortNLMrot, 
    #                      areaAvg_negVortNLMrot],axis=0)
    
    return returnArr

if __name__ == '__main__':
    #rootFolder = '/discover/nobackup/projects/mesocean/srai/runGeos7dayAvg/'

    
    #ds_masks = Dataset(rootFolder + 'tavgInput/RegionMasks.nc')
    #GridFile = rootFolder + 'tavgInput/newGlobalGrid_tripolePOP_0.1deg.nc'
    
    rootFolder = '/pscratch/sd/s/srai/nccs/runGeos7dayAvg/'
    ds_masks = Dataset(rootFolder + 'test_7dayAvg_NanLand_test/input/RegionMasks.nc')
    GridFile = rootFolder + 'test_7dayAvg_NanLand_test/input/tripoleGridCreated.nc'
    kmtSatFile = 'test_7dayAvg_NanLand_test/input/kmtSatToCESM.nc'
    keAndEkeFile = rootFolder + 'test_7dayAvg_NanLand_test/input/KEandEKE_pop.nc'
    

    kmtDS = Dataset(rootFolder + kmtSatFile)
    KMT = 1-np.roll(np.array(kmtDS.variables['KMTsat']), 1800, axis=1)
    
    keDS = Dataset(keAndEkeFile)
    KE = np.array(keDS.variables['KE'][:,:])/10000
    EKE = np.array(keDS.variables['EKE'][:,:])/10000
    highKE_mask = (KE + EKE) > 0.04    
    
    gridDS = Dataset(GridFile)
    UAREA = np.array(gridDS.variables['UAREA'])
    ULAT = np.array(gridDS.variables['ULAT'])
    ULONG = np.array(gridDS.variables['ULONG'])
    #KMT = np.array(gridDS.variables['KMT'])
    ylen, xlen = np.shape(UAREA)
    DX = np.array(gridDS.variables['DXU'])
    DY = np.array(gridDS.variables['DYU'])

    mask_20_40N = np.logical_and(ULAT > 20, ULAT < 40)
    mask_5S_5N = np.logical_and(ULAT > -5, ULAT < 5)
    mask_15_40S = np.logical_and(ULAT > -40, ULAT < -15)
    landMask = KMT < 1
    #highKE_mask = np.array(ds_masks['HighEKE'][:, :])

    landMaskWithHiKE = np.logical_or(landMask, highKE_mask)
    UAREA[landMaskWithHiKE] = float('nan')
    #plt.pcolormesh(landMaskWithHiKE)
    #plt.colorbar()
    #plt.show()


    # ellList = np.arange(60, 450, 20)
    # ellList = np.append(np.array([50]), ellList, axis=0)
    # ellList = np.append(ellList, np.arange(500, 1050, 100), axis=0)
    ellList = np.array([100])
    nell = len(ellList)
    ndays = 365

    for ellIDX in range(nell):
        ell = ellList[ellIDX]
        ellFold = rootFolder + 'test_7dayAvg_NanLand_test/Output/{0:d}km_smth/'.format(ell)

        writeFileName = ellFold + \
            'NLMrotAndNLMstrLatBandSeries_{0:04d}km.nc'.format(
                ell)

        fileName1 = 'NLMonPosAndNegVort_{0:04d}km.nc'.format(
            ell)
        fileName2 = 'NLMon45degAndMinus45degStrain_{0:04d}km.nc'.format(
            ell)

        print('working in file {0:s} and {1:s} for ell {2:d}km'.format(
            fileName1, fileName2, ell))

        ds1 = Dataset(ellFold + fileName1)
        ds2 = Dataset(ellFold + fileName2)

        NLMstrNrot_20_40N = np.zeros((ndays,4), dtype=float)
        NLMstrNrot_15_40S = np.zeros((ndays,4), dtype=float)
        NLMstrNrot_5S_5N = np.zeros((ndays,4), dtype=float)

        for day in range(ndays):
            print('reading day', day)
            posVortNLMrot = np.array(ds1.variables['posVortNLMrot'][day, :, :])
            negVortNLMrot = np.array(ds1.variables['negVortNLMrot'][day, :, :])
            posThetaNLMstr = np.array(ds2.variables['posThetaNLMstr'][day, :, :])
            negThetaNLMstr = np.array(ds2.variables['negThetaNLMstr'][day, :, :])

            posVortNLMrot[landMaskWithHiKE] = float('nan')
            negVortNLMrot[landMaskWithHiKE] = float('nan')
            posThetaNLMstr[landMaskWithHiKE] = float('nan')
            negThetaNLMstr[landMaskWithHiKE] = float('nan')

            mask = mask_20_40N.copy()
            NLMstrNrot_20_40N[day,:] = getAreaAvgMaskedVal(posThetaNLMstr,
                                                    negThetaNLMstr,
                                                    posVortNLMrot,
                                                    negVortNLMrot,
                                                    UAREA, mask)
            
            mask = mask_15_40S.copy()
            NLMstrNrot_15_40S[day,:] = getAreaAvgMaskedVal(posThetaNLMstr,
                                                    negThetaNLMstr,
                                                    posVortNLMrot,
                                                    negVortNLMrot,
                                                    UAREA, mask)
            
            mask = mask_5S_5N.copy()
            NLMstrNrot_5S_5N[day,:] = getAreaAvgMaskedVal(posThetaNLMstr,
                                                    negThetaNLMstr,
                                                    posVortNLMrot,
                                                    negVortNLMrot,
                                                    UAREA, mask)
        
        
        wds = Dataset(writeFileName, 'w', format='NETCDF4')

        wds.createDimension('time', ndays)

        createVariableAndWrite(wds, 'posThetaNLMstr_20_40N', 'ergs/cm^2/sec',
                               'NLM_str on thetaVel > 0',
                               ('time'),
                               float,
                               NLMstrNrot_20_40N[:,0])
        
        createVariableAndWrite(wds, 'negThetaNLMstr_20_40N', 'ergs/cm^2/sec',
                               'NLM_str on thetaVel < 0',
                               ('time'),
                               float,
                               NLMstrNrot_20_40N[:, 1])
        
        createVariableAndWrite(wds, 'posVortNLMrot_20_40N', 'ergs/cm^2/sec',
                               'NLM_rot on vorticity > 0',
                               ('time'),
                               float,
                               NLMstrNrot_20_40N[:, 2])

        createVariableAndWrite(wds, 'negVortNLMrot_20_40N', 'ergs/cm^2/sec',
                               'NLM_rot on vorticity < 0',
                               ('time'),
                               float,
                               NLMstrNrot_20_40N[:, 3])
        
        createVariableAndWrite(wds, 'posThetaNLMstr_5S_5N', 'ergs/cm^2/sec',
                               'NLM_str on thetaVel > 0',
                               ('time'),
                               float,
                               NLMstrNrot_5S_5N[:, 0])

        createVariableAndWrite(wds, 'negThetaNLMstr_5S_5N', 'ergs/cm^2/sec',
                               'NLM_str on thetaVel < 0',
                               ('time'),
                               float,
                               NLMstrNrot_5S_5N[:, 1])

        createVariableAndWrite(wds, 'posVortNLMrot_5S_5N', 'ergs/cm^2/sec',
                               'NLM_rot on vorticity > 0',
                               ('time'),
                               float,
                               NLMstrNrot_5S_5N[:, 2])

        createVariableAndWrite(wds, 'negVortNLMrot_5S_5N', 'ergs/cm^2/sec',
                               'NLM_rot on vorticity < 0',
                               ('time'),
                               float,
                               NLMstrNrot_5S_5N[:, 3])

        createVariableAndWrite(wds, 'posThetaNLMstr_15_40S', 'ergs/cm^2/sec',
                               'NLM_str on thetaVel > 0',
                               ('time'),
                               float,
                               NLMstrNrot_15_40S[:, 0])

        createVariableAndWrite(wds, 'negThetaNLMstr_15_40S', 'ergs/cm^2/sec',
                               'NLM_str on thetaVel < 0',
                               ('time'),
                               float,
                               NLMstrNrot_15_40S[:, 1])

        createVariableAndWrite(wds, 'posVortNLMrot_15_40S', 'ergs/cm^2/sec',
                               'NLM_rot on vorticity > 0',
                               ('time'),
                               float,
                               NLMstrNrot_15_40S[:, 2])

        createVariableAndWrite(wds, 'negVortNLMrot_15_40S', 'ergs/cm^2/sec',
                               'NLM_rot on vorticity < 0',
                               ('time'),
                               float,
                               NLMstrNrot_15_40S[:, 3])

        wds.close()
            
            


