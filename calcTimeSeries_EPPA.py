from netCDF4 import Dataset
import numpy as np

fldLoc = '/pscratch/sd/s/srai/nccs/runGeos7dayAvg/test_7dayAvg_NanLand_test/input/'
gridFile = fldLoc + 'tripoleGridCreated.nc'
maskFile = fldLoc + 'RegionMasks.nc'
satMaskFile = fldLoc + 'kmtSatToCESM.nc'

gridDS2 = Dataset(satMaskFile)
satKMT = np.array(gridDS2.variables['KMTsat'])
satLandMask = np.roll(satKMT > 0,-1800, axis=1)

# globalAreaDS = Dataset(fldLoc + '/GlbAreaList.nc')
# glbArea = np.array(globalAreaDS.variables['newGlobalGrid_tripolePOP_0.1deg.nc'])

# ellList = np.arange(60, 450, 20)
# ellList = np.append(np.array([50]), ellList, axis=0)
# ellList = np.append(ellList, np.arange(500, 1050, 100), axis=0)

ellList = [100]

nell = len(ellList)
ndays = 365
ylen = 2400
xlen = 3600

gridDS = Dataset(gridFile)
maskDS = Dataset(maskFile)

regionID = ['Gulf',
            'Kuroshio',
            'BrazilACC',
            'Agulhas',
            'SouthernOcean',
            'ACC',
            'Equator',
            'Global',
            'HighEKE']

nregions = len(regionID)

masks = np.ones((nregions, ylen, xlen), dtype=float)
for i in range(nregions-1):
    masks[i, :, :] = np.array(maskDS.variables[regionID[i]][:, :])


TPPA = np.zeros((nell, nregions, ndays), dtype=float)
MPPA = np.zeros((nell, nregions, ndays), dtype=float)
EPPA = np.zeros((nell, nregions, ndays), dtype=float)

TPPA_PerArea = np.zeros((nell, nregions, ndays), dtype=float)
MPPA_PerArea = np.zeros((nell, nregions, ndays), dtype=float)
EPPA_PerArea = np.zeros((nell, nregions, ndays), dtype=float)

UAREA = np.array(gridDS.variables['UAREA'])
UAREA[satLandMask] = 0.0
KMT = np.array(gridDS.variables['KMT'])
landMask = KMT < 1

for ellIDX in range(nell):
    ell = ellList[ellIDX]
    print('ell = {0:d}'.format(ell))
    fileName = '{0:d}km_smth/filtered_{0:04d}.nc'.format(ell)
    
    ds = Dataset(fileName)

    for dumIDX in range(ndays):
        print('day', dumIDX)
        this_EPPA = np.array(
            ds.variables['EddyPowerPerArea'][dumIDX, :, :])
        this_MPPA = np.array(
            ds.variables['MeanPowerPerArea'][dumIDX, :, :])
        this_TPPA = np.array(
            ds.variables['TotalPowerPerArea'][dumIDX, :, :])
        
        this_EPPA[abs(this_EPPA)>1e10] = 0.0
        this_MPPA[abs(this_MPPA)>1e10] = 0.0
        this_TPPA[abs(this_TPPA)>1e10] = 0.0

        this_EPPA[satLandMask] = 0.0
        this_MPPA[satLandMask] = 0.0
        this_TPPA[satLandMask] = 0.0
        
        for regionIDX in range(nregions):
            mskUAREA = masks[regionIDX, :, :] * UAREA
            totUAREA = np.sum(mskUAREA)

            if regionID[regionIDX] == 'Global':
                mask = np.isnan(this_EPPA) + landMask
                mskUAREA = np.ma.array(UAREA, mask=mask, fill_value=0).filled()
                totUAREA = np.nansum(mskUAREA)

            TPPA[ellIDX, regionIDX, dumIDX] = np.nansum(
                mskUAREA[:, :] * this_TPPA[:, :])

            MPPA[ellIDX, regionIDX, dumIDX] = np.nansum(
                mskUAREA[:, :] * this_MPPA[:, :])

            EPPA[ellIDX, regionIDX, dumIDX] = np.nansum(
                mskUAREA[:, :] * this_EPPA[:, :])

            TPPA_PerArea[ellIDX, regionIDX,
                            dumIDX] = TPPA[ellIDX, regionIDX, dumIDX]/totUAREA

            MPPA_PerArea[ellIDX, regionIDX,
                            dumIDX] = MPPA[ellIDX, regionIDX, dumIDX]/totUAREA

            EPPA_PerArea[ellIDX, regionIDX,
                            dumIDX] = EPPA[ellIDX, regionIDX, dumIDX]/totUAREA


outfileName = 'allEllallRegionsEP.nc'

wds = Dataset(outfileName, 'w', format='NETCDF4')

wds.createDimension('ell', nell)
wds.createDimension('region', nregions)
wds.createDimension('time', None)

wds_ell = wds.createVariable('ell', int, ('ell'))
wds_ell.units = 'km'
wds_ell.long_name = 'filterlength'
wds_ell[:] = ellList[:]

wds_regions = wds.createVariable('region', str, ('region'))
# wds_regions.units = ''
wds_regions.long_name = 'region masked'
wds_regions[:] = np.array(regionID, dtype='object')

wds_time = wds.createVariable('time', float, ('time'))
wds_time.units = 'days since 0000-1-1 00:00 UTC'
wds_time.long_name = 'date'
wds_time[:] = np.arange(365) + 365*49 + 3

wdsVar_TPPA = wds.createVariable(
    'TPPA', float, ('ell', 'region', 'time'))
wdsVar_TPPA.units = 'ergs/(m^-2 sec)'
wdsVar_TPPA.long_name = 'TP^(Cg) Per Area'
wdsVar_TPPA[:, :, :] = TPPA_PerArea[:, :, :]

wdsVar_MPPA = wds.createVariable(
    'MPPA', float, ('ell', 'region', 'time'))
wdsVar_MPPA.units = 'ergs/(m^-2 sec)'
wdsVar_MPPA.long_name = 'MP^(Cg) Per Area'
wdsVar_MPPA[:, :, :] = MPPA_PerArea[:, :, :]

wdsVar_EPPA = wds.createVariable(
    'EPPA', float, ('ell', 'region', 'time'))
wdsVar_EPPA.units = 'ergs/(m^-2 sec)'
wdsVar_EPPA.long_name = 'EP^(Cg) Per Area'
wdsVar_EPPA[:, :, :] = EPPA_PerArea[:, :, :]


# wdsVar_TPPA = wds.createVariable(
#     'TPPA', float, ('ell', 'region', 'time'))
# wdsVar_TPPA.units = 'watts'
# wdsVar_TPPA.long_name = 'total modelled EP^(Cg)'
# wdsVar_TPPA[:, :, :] = TPPA[:, :, :]

# wdsVar_NLM2Order_tot = wds.createVariable(
#     'NLM2Order_tot', float, ('ell', 'region', 'time'))
# wdsVar_NLM2Order_tot.units = 'watts'
# wdsVar_NLM2Order_tot.long_name = 'total modelled EP^(Cg)'
# wdsVar_NLM2Order_tot[:, :, :] = NLM2Order_tot[:, :, :]

# wdsVar_MPPA = wds.createVariable(
#     'MPPA', float, ('ell', 'region', 'time'))
# wdsVar_MPPA.units = 'watts'
# wdsVar_MPPA.long_name = 'vortical component of modelled EP^(Cg)'
# wdsVar_MPPA[:, :, :] = MPPA[:, :, :]

# wdsVar_EPPA = wds.createVariable(
#     'EPPA', float, ('ell', 'region', 'time'))
# wdsVar_EPPA.units = 'watts'
# wdsVar_EPPA.long_name = 'strain component of modelled EP^(Cg)'
# wdsVar_EPPA[:, :, :] = EPPA[:, :, :]


wds.close()

