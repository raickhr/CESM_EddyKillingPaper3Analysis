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


NLM_tot = np.zeros((nell, nregions, ndays), dtype=float)
NLM2Order_tot = np.zeros((nell, nregions, ndays), dtype=float)
NLM_rot = np.zeros((nell, nregions, ndays), dtype=float)
NLM_str = np.zeros((nell, nregions, ndays), dtype=float)

NLM_tot_PerArea = np.zeros((nell, nregions, ndays), dtype=float)
NLM2Order_tot_PerArea = np.zeros((nell, nregions, ndays), dtype=float)
NLM_rot_PerArea = np.zeros((nell, nregions, ndays), dtype=float)
NLM_str_PerArea = np.zeros((nell, nregions, ndays), dtype=float)

UAREA = np.array(gridDS.variables['UAREA'])
satLandMask[satLandMask] = 0.0
KMT = np.array(gridDS.variables['KMT'])
landMask = KMT < 1

for ellIDX in range(nell):
    ell = ellList[ellIDX]
    print('ell = {0:d}'.format(ell))
    fileName = '{0:d}km_smth/NLmodelEP_{0:04d}km.nc'.format(ell)
    #slopeFile = '{0:d}km/slopeAndCorr2D_{0:04d}km.nc'.format(ell)

    ds = Dataset(fileName)
    #slopeDS = Dataset(slopeFile)

    slope = 1 #np.array(slopeDS.variables['slope'])
    
    for dumIDX in range(ndays):
        this_NLM_tot = np.array(
            ds.variables['NLmodel_EPCg'][dumIDX, :, :])/slope
        this_NLM2Order_tot = np.array(
            ds.variables['NLmodel2_EPCg'][dumIDX, :, :])/slope
        this_NLM_rot = np.array(
            ds.variables['NLmodel_EPCg_rot'][dumIDX, :, :])/slope
        this_NLM_str = np.array(
            ds.variables['NLmodel_EPCg_strain'][dumIDX, :, :])/slope
        
        for regionIDX in range(nregions):
            mskUAREA = masks[regionIDX, :, :] * UAREA
            totUAREA = np.sum(mskUAREA)

            if regionID[regionIDX] == 'Global':
                mask = np.isnan(this_NLM2Order_tot) + landMask
                mskUAREA = np.ma.array(UAREA,mask=mask, fill_value=0).filled()
                totUAREA = np.nansum(mskUAREA)

            NLM_tot[ellIDX, regionIDX, dumIDX] = np.nansum(
                mskUAREA[:, :] * this_NLM_tot[:, :])
            
            NLM2Order_tot[ellIDX, regionIDX, dumIDX] = np.nansum(
                mskUAREA[:, :] * this_NLM2Order_tot[:, :])
            
            NLM_rot[ellIDX, regionIDX, dumIDX] = np.nansum(
                mskUAREA[:, :] * this_NLM_rot[:, :])
            
            NLM_str[ellIDX, regionIDX, dumIDX] = np.nansum(
                mskUAREA[:, :] * this_NLM_str[:, :])

            NLM_tot_PerArea[ellIDX, regionIDX,
                            dumIDX] = NLM_tot[ellIDX, regionIDX, dumIDX]/totUAREA
            
            NLM2Order_tot_PerArea[ellIDX, regionIDX,
                                    dumIDX] = NLM2Order_tot[ellIDX, regionIDX, dumIDX]/totUAREA
            
            NLM_rot_PerArea[ellIDX, regionIDX,
                            dumIDX] = NLM_rot[ellIDX, regionIDX, dumIDX]/totUAREA
            
            NLM_str_PerArea[ellIDX, regionIDX,
                            dumIDX] = NLM_str[ellIDX, regionIDX, dumIDX]/totUAREA


outfileName = 'allEllallRegionsEP_ModelRotAndStrComp.nc'

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

wdsVar_NLM_tot = wds.createVariable(
    'NLM_tot_PerArea', float, ('ell', 'region', 'time'))
wdsVar_NLM_tot.units = 'ergs/(m^-2 sec)'
wdsVar_NLM_tot.long_name = 'total modelled EP^(Cg) Per Area'
wdsVar_NLM_tot[:, :, :] = NLM_tot_PerArea[:, :, :]

wdsVar_NLM2Order_tot = wds.createVariable(
    'NLM2Order_tot_PerArea', float, ('ell', 'region', 'time'))
wdsVar_NLM2Order_tot.units = 'ergs/(m^-2 sec)'
wdsVar_NLM2Order_tot.long_name = 'total modelled EP^(Cg) Per Area'
wdsVar_NLM2Order_tot[:, :, :] = NLM2Order_tot_PerArea[:, :, :]

wdsVar_NLM_rot = wds.createVariable(
    'NLM_rot_PerArea', float, ('ell', 'region', 'time'))
wdsVar_NLM_rot.units = 'ergs/(m^-2 sec)'
wdsVar_NLM_rot.long_name = 'vortical component of modelled EP^(Cg) Per Area'
wdsVar_NLM_rot[:, :, :] = NLM_rot_PerArea[:, :, :]

wdsVar_NLM_str = wds.createVariable(
    'NLM_str_PerArea', float, ('ell', 'region', 'time'))
wdsVar_NLM_str.units = 'ergs/(m^-2 sec)'
wdsVar_NLM_str.long_name = 'strain component of modelled EP^(Cg) Per Area'
wdsVar_NLM_str[:, :, :] = NLM_str_PerArea[:, :, :]


# wdsVar_NLM_tot = wds.createVariable(
#     'NLM_tot', float, ('ell', 'region', 'time'))
# wdsVar_NLM_tot.units = 'watts'
# wdsVar_NLM_tot.long_name = 'total modelled EP^(Cg)'
# wdsVar_NLM_tot[:, :, :] = NLM_tot[:, :, :]

# wdsVar_NLM2Order_tot = wds.createVariable(
#     'NLM2Order_tot', float, ('ell', 'region', 'time'))
# wdsVar_NLM2Order_tot.units = 'watts'
# wdsVar_NLM2Order_tot.long_name = 'total modelled EP^(Cg)'
# wdsVar_NLM2Order_tot[:, :, :] = NLM2Order_tot[:, :, :]

# wdsVar_NLM_rot = wds.createVariable(
#     'NLM_rot', float, ('ell', 'region', 'time'))
# wdsVar_NLM_rot.units = 'watts'
# wdsVar_NLM_rot.long_name = 'vortical component of modelled EP^(Cg)'
# wdsVar_NLM_rot[:, :, :] = NLM_rot[:, :, :]

# wdsVar_NLM_str = wds.createVariable(
#     'NLM_str', float, ('ell', 'region', 'time'))
# wdsVar_NLM_str.units = 'watts'
# wdsVar_NLM_str.long_name = 'strain component of modelled EP^(Cg)'
# wdsVar_NLM_str[:, :, :] = NLM_str[:, :, :]


wds.close()

