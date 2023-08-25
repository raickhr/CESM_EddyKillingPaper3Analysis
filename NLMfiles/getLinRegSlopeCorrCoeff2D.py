from netCDF4 import Dataset
import numpy as np
from scipy import stats

n = 365
Ylen = 2400
Xlen = 3600
ell = 100

GridFile = '/discover/nobackup/projects/mesocean/srai/runGeos7dayAvg/test_7dayAvg_NanLand_test/input/tripoleGridCreated.nc'
ds = Dataset(GridFile)
KMT = np.array(ds.variables['KMT'])
landMask = KMT < 1


sum_y = np.zeros((Ylen, Xlen), dtype=float)
sum_x = np.zeros((Ylen, Xlen), dtype=float)
sum_xy = np.zeros((Ylen, Xlen), dtype=float)
sum_xsq = np.zeros((Ylen, Xlen), dtype=float)
sum_ysq = np.zeros((Ylen, Xlen), dtype=float)

fldLoc = '/discover/nobackup/projects/mesocean/srai/runGeos7dayAvg/test_7dayAvg_NanLand_test/Output/100km_smth/'

print('calculating mean values')
for t in range(n):
    print('at file {0:d}'.format(t+1))
    yfileName = '{0:03d}.nc_{1:04d}_Filtered.nc'.format(t+1, ell)
    xfileName = '{0:03d}_NLmodelEP_{1:04d}km.nc'.format(t+1, ell)
    xds = Dataset(fldLoc + xfileName)
    yds = Dataset(fldLoc + yfileName)

    y = np.array(yds.variables['EddyPowerPerArea'][0,:,:])
    y[landMask] = float('nan')

    x = np.array(xds.variables['NLmodel_EPCg'][0,:,:]) + \
        np.array(xds.variables['NLmodel2_EPCg'][0,:,:])
    x[landMask] = float('nan')

    sum_x += x
    sum_y += y
    sum_xy += x*y
    sum_xsq += x**2
    sum_ysq += y**2

##################### mean #####################

mean_x = sum_x/n
mean_y = sum_y/n

mean_xy = sum_xy/n
mean_xsq = sum_xsq/n
mean_ysq = sum_ysq/n

##################### std dev #################

var_x = mean_xsq - mean_x**2
var_y = mean_ysq - mean_y**2

sigma_x = np.sqrt(n/(n-1)*var_x)
sigma_y = np.sqrt(n/(n-1)*var_y)

###################### intercept ################

nume = sum_y*sum_xsq - sum_x*sum_xy
deno = n * sum_xsq - sum_x**2
c = nume/deno

###################### slope ################

nume = n*sum_xy - (sum_x * sum_y)
deno = n * sum_xsq - sum_x**2
m = nume/deno


###################### correlation coeff ################

nume = ((mean_xy) - (mean_x)*(mean_y))
deno = np.sqrt((mean_xsq - (mean_x)**2)*(mean_ysq - (mean_y)**2))
rxy = nume/deno

###################### standard err of slope ################


sum_yMinusYhat_sq = np.zeros((Ylen, Xlen), dtype=float)

print('\n\ncalculating std error for slope')
for t in range(n):
    print('at file {0:d}'.format(t+1))
    yfileName = '{0:03d}.nc_{1:04d}_Filtered.nc'.format(t+1, ell)
    xfileName = '{0:03d}_NLmodelEP_{1:04d}km.nc'.format(t+1, ell)
    xds = Dataset(fldLoc + xfileName)
    yds = Dataset(fldLoc + yfileName)

    y = np.array(yds.variables['EddyPowerPerArea'][0, :, :])
    y[landMask] = float('nan')

    x = np.array(xds.variables['NLmodel_EPCg'][0, :, :]) + \
        np.array(xds.variables['NLmodel2_EPCg'][0, :, :])
    x[landMask] = float('nan')

    Yhat = c + m * x
    yMinusYhat_sq = (y - Yhat)**2

    sum_yMinusYhat_sq += yMinusYhat_sq

nume = sum_yMinusYhat_sq
deno = (n-2) * (mean_xsq - mean_x**2)
std_err = np.sqrt(nume/deno)

prob = 0.95

t_star = stats.t.cdf(1 - (1-prob)/2, n - 1)
std_err *= t_star * sigma_x/sigma_y



######################  writing netcdf file ######################

writeFileName = '{0:d}km_smth/'.format(ell) + 'slopeAndCorr2D_{0:04d}km.nc'.format(ell)

wds = Dataset(writeFileName, 'w', 'NETCDF4')
wds.createDimension('Y', Ylen)
wds.createDimension('X', Xlen)

wds_slope = wds.createVariable('slope', float, ('Y', 'X'))
wds_intercept = wds.createVariable('intercept', float, ('Y', 'X'))
wds_sigma_x = wds.createVariable('sigma_x', float, ('Y', 'X'))
wds_sigma_y = wds.createVariable('sigma_y', float, ('Y', 'X'))
wds_std_err = wds.createVariable('std_err', float, ('Y', 'X'))
wds_corr_coeff = wds.createVariable('corr_coeff', float, ('Y', 'X'))

wds_slope[:, :] = m[:,:]
wds_intercept[:, :] = c[:,:]
wds_sigma_x[:,:] = sigma_x[:,:]
wds_sigma_y[:,:] = sigma_y[:,:]
wds_std_err[:, :] = std_err[:,:]
wds_corr_coeff[:, :] = rxy[:,:]

wds.close()





