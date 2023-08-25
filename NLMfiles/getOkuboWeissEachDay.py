from netCDF4 import Dataset
import numpy as np

def calcGradient(Array2D, DXU, DYU):
    p = -1
    m = 1
    x = 1
    y = 0
    topRight = 0.25*(Array2D + np.roll(Array2D, p, axis=x) +
                     np.roll(Array2D, (p, p), axis=(y, x)) +
                     np.roll(Array2D, p, axis=y))

    topLeft = 0.25*(Array2D + np.roll(Array2D, p, axis=y) +
                    np.roll(Array2D, (p, m), axis=(y, x)) +
                    np.roll(Array2D, m, axis=x))

    bottomRight = 0.25*(Array2D + np.roll(Array2D, p, axis=x) +
                        np.roll(Array2D, (m, p), axis=(y, x)) +
                        np.roll(Array2D, m, axis=y))

    bottomLeft = 0.25*(Array2D + np.roll(Array2D, m, axis=x) +
                       np.roll(Array2D, (m, m), axis=(y, x)) +
                       np.roll(Array2D, m, axis=y))

    gradx = 0.5*(topRight + bottomRight - topLeft - bottomLeft)/DXU
    grady = 0.5*(topRight + topLeft - bottomRight - bottomLeft)/DYU

    return gradx, grady


def getCurlZ(u, v, DX, DY):
    dx_u, dy_u = calcGradient(u, DX, DY)
    dx_v, dy_v = calcGradient(v, DX, DY)
    zcurl = dx_v - dy_u
    return zcurl


def getStrainTensor(u, v, DX, DY):
    dx_u, dy_u = calcGradient(u, DX, DY)
    dx_v, dy_v = calcGradient(v, DX, DY)
    S11 = dx_u
    S12 = 0.5 * (dx_v + dy_u)
    S21 = S12
    S22 = dy_v
    return S11, S12, S21, S22


def getOkuboWeissParameter(u, v, DX, DY):
    dx_u, dy_u = calcGradient(u, DX, DY)
    dx_v, dy_v = calcGradient(v, DX, DY)
    normalStrain = dx_u - dy_v
    shearStrain = dx_v + dy_u
    omega = dx_v - dy_u
    W = normalStrain**2 + shearStrain**2 - omega**2
    return W


def rangePlusMinus90(theta):
    theta[theta > 90] -= 180
    theta[theta < -90] += 180
    return theta


def getStrainEigValEigVec(S11, S12, S21, S22):
    strainTensor = np.array([[S11, S12],
                             [S21, S22]], dtype=float)

    w, v = np.linalg.eig(strainTensor.transpose())

    y1 = v[:, :, 0, 1].transpose()
    x1 = v[:, :, 0, 0].transpose()

    y2 = v[:, :, 1, 1].transpose()
    x2 = v[:, :, 1, 0].transpose()

    theta1 = np.degrees(np.arctan2(y1, x1))
    theta2 = np.degrees(np.arctan2(y2, x2))

    eig1 = w[:, :, 0].transpose()
    eig2 = w[:, :, 1].transpose()

    mask = eig2 > eig1

    dum1 = eig1.copy()
    dum2 = eig2.copy()

    eig1[mask] = dum2[mask]
    eig2[mask] = dum1[mask]

    del dum1, dum2

    dum1 = theta1.copy()
    dum2 = theta2.copy()

    theta1[mask] = dum2[mask]
    theta2[mask] = dum1[mask]

    del dum1, dum2

    return eig1, eig2, -rangePlusMinus90(theta1), -rangePlusMinus90(theta2)


def getTransMat(S11, S12, S21, S22,
                T11, T12, T21, T22):
    strainTensor = np.array([[S11, S12],
                             [S21, S22]], dtype=float)

    TauTensor = np.array([[T11, T12],
                          [T21, T22]], dtype=float)

    # invStrainTensor = np.linalg.inv(strainTensor.transpose())

    invTauTensor = np.linalg.inv(TauTensor.transpose())

    transFormMat = np.matmul(strainTensor.transpose(), invTauTensor)

    # transFormMat = np.matmul(TauTensor.transpose(), invStrainTensor)

    return transFormMat.transpose()


def createVariableAndWrite(wds, varName, varUnits, varLongName, varDimension, varType, varArr):
    cdfVar = wds.createVariable(varName, varType, varDimension)
    cdfVar.units = varUnits
    cdfVar.long_name = varLongName
    cdfVar[:] = varArr[:]


if __name__ == '__main__':
    #rootFolder = '/discover/nobackup/projects/mesocean/srai/runGeos7dayAvg/'

    #ds_masks = Dataset(rootFolder + 'tavgInput/RegionMasks.nc')
    #GridFile = rootFolder + 'tavgInput/newGlobalGrid_tripolePOP_0.1deg.nc'
    
    rootFolder = '/pscratch/sd/s/srai/nccs/runGeos7dayAvg/test_7dayAvg_NanLand_test/'
    ds_masks = Dataset(rootFolder + 'input/RegionMasks.nc')
    GridFile = rootFolder + 'input/tripoleGridCreated.nc'

    gridDS = Dataset(GridFile)
    UAREA = np.array(gridDS.variables['UAREA'])
    KMT = np.array(gridDS.variables['KMT'])
    ylen, xlen = np.shape(UAREA)
    DX = np.array(gridDS.variables['DXU'])
    DY = np.array(gridDS.variables['DYU'])
    DX[DX == 0] = 7e8  ## approx. raidus of earth
    DY[DY == 0] = 7e8  ## approx. raidus of earth

    # ellList = np.arange(60, 450, 20)
    # ellList = np.append(np.array([50]), ellList, axis=0)
    # ellList = np.append(ellList, np.arange(500, 1050, 100), axis=0)
    ellList = np.array([100])
    nell = len(ellList)
    ndays = 365

    landMask = KMT < 1
    
    for ellIDX in range(nell):
        ell = ellList[ellIDX]
        print('ell = ', ell)
        ellFold = rootFolder + 'Output/{0:d}km_smth/'.format(ell)


        fileName = 'filtered_{0:04}.nc'.format(ell)
        ds = Dataset(ellFold + fileName)

        TIME = ds.variables['time']

        print('all variables read')

        for timeIdx in range(326,len(TIME)):
            writeFileName = ellFold + \
                'okuboWeissAndStrainDirec_{0:04d}km_Day_{1:03d}.nc'.format(
                    ell,timeIdx+1)
            wds = Dataset(writeFileName, 'w', format='NETCDF4')

            wds.createDimension('Y', ylen)
            wds.createDimension('X', xlen)
            wds.createDimension('time', None)

            cdftime = wds.createVariable('time', float, ('time'))
            cdftime.units = TIME.units
            timeArr = np.array(TIME)
            cdftime[:] = timeArr[timeIdx]

            okuboWeiss_vel = np.zeros((ylen, xlen), dtype=float)
            eig1_tau = np.zeros((ylen, xlen), dtype=float)
            eig2_tau = np.zeros((ylen, xlen), dtype=float)
            theta1_tau = np.zeros((ylen, xlen), dtype=float)
            theta2_tau = np.zeros((ylen, xlen), dtype=float)
            eig1_vel = np.zeros((ylen, xlen), dtype=float)
            eig2_vel = np.zeros((ylen, xlen), dtype=float)
            theta1_vel = np.zeros((ylen, xlen), dtype=float)
            theta2_vel = np.zeros((ylen, xlen), dtype=float)
            omega_tau = np.zeros((ylen, xlen), dtype=float)
            omega_vel = np.zeros((ylen, xlen), dtype=float)
            omega_tau = np.zeros((ylen, xlen), dtype=float)
            omega_vel = np.zeros((ylen, xlen), dtype=float)

            TAUX = np.array(ds.variables['SUF'][timeIdx])
            TAUY = np.array(ds.variables['SVF'][timeIdx])
            UVEL = np.array(ds.variables['UVEL'][timeIdx])
            VVEL = np.array(ds.variables['VVEL'][timeIdx])

            TAUX[abs(TAUX)>1e5] = float('nan')
            TAUY[abs(TAUX)>1e5] = float('nan')
            UVEL[abs(TAUX)>1e5] = float('nan')
            VVEL[abs(TAUX)>1e5] = float('nan')

            print('in day {0:d}'.format(timeIdx + 1))
            okuboWeiss_vel[:, :] = getOkuboWeissParameter(
                UVEL[:, :], VVEL[:, :], DX, DY)
            
            S11_tau, S12_tau, S21_tau, S22_tau = getStrainTensor(
                TAUX[:, :], TAUY[:, :], DX, DY)

            S11_tau[landMask] = 0.0
            S12_tau[landMask] = 0.0
            S21_tau[landMask] = 0.0
            S22_tau[landMask] = 0.0

            S11_tau[np.isinf(S11_tau)] = 0.0
            S12_tau[np.isinf(S12_tau)] = 0.0
            S21_tau[np.isinf(S21_tau)] = 0.0
            S22_tau[np.isinf(S22_tau)] = 0.0

            S11_tau[np.isnan(S11_tau)] = 0.0
            S12_tau[np.isnan(S12_tau)] = 0.0
            S21_tau[np.isnan(S21_tau)] = 0.0
            S22_tau[np.isnan(S22_tau)] = 0.0

            eig1_tau[:, :], eig2_tau[:, :], theta1_tau[:, :], theta2_tau[:, :] = getStrainEigValEigVec(
                S11_tau, S12_tau, S21_tau, S22_tau)
            del S11_tau, S12_tau, S21_tau, S22_tau

            eig1_tau[landMask] = float('nan')
            eig2_tau[landMask] = float('nan')
            theta1_tau[landMask] = float('nan')
            theta2_tau[landMask] = float('nan')

            S11_vel, S12_vel, S21_vel, S22_vel = getStrainTensor(
                UVEL[:, :], VVEL[:, :], DX, DY)

            S11_vel[landMask] = 0.0
            S12_vel[landMask] = 0.0
            S21_vel[landMask] = 0.0
            S22_vel[landMask] = 0.0

            S11_vel[np.isinf(S11_vel)] = 0.0
            S12_vel[np.isinf(S12_vel)] = 0.0
            S21_vel[np.isinf(S21_vel)] = 0.0
            S22_vel[np.isinf(S22_vel)] = 0.0

            S11_vel[np.isnan(S11_vel)] = 0.0
            S12_vel[np.isnan(S12_vel)] = 0.0
            S21_vel[np.isnan(S21_vel)] = 0.0
            S22_vel[np.isnan(S22_vel)] = 0.0

            eig1_vel[:, :], eig2_vel[:, :], theta1_vel[:, :], theta2_vel[:, :] = getStrainEigValEigVec(
                S11_vel, S12_vel, S21_vel, S22_vel)

            del S11_vel, S12_vel, S21_vel, S22_vel

            eig1_vel[landMask] = float('nan')
            eig2_vel[landMask] = float('nan')
            theta1_vel[landMask] = float('nan')
            theta2_vel[landMask] = float('nan')

            omega_tau[:, :] = getCurlZ(
                TAUX[:, :], TAUY[:, :], DX, DY)
            omega_vel[:, :] = getCurlZ(
                UVEL[:, :], VVEL[:, :], DX, DY)
            omega_tau[landMask] = float('nan')
            omega_vel[landMask] = float('nan')

            createVariableAndWrite(wds, 'okuboWeiss_vel', 'N/A',
                                'okubo weiss parameter for ocean current',
                                ('time', 'Y', 'X'),
                                float,
                                np.array([okuboWeiss_vel],dtype=float))

            createVariableAndWrite(wds, 'eig1_tau', 'dynes/cm^3',
                                'first eigenvalue of the symmetric part of stress gradient tensor',
                                ('time', 'Y', 'X'),
                                float,
                                np.array([eig1_tau],dtype=float))

            createVariableAndWrite(wds, 'eig2_tau', 'dynes/cm^3',
                                'second eigenvalue of the symmetric part of stress gradient tensor',
                                ('time', 'Y', 'X'),
                                float,
                                np.array([eig2_tau],dtype=float))

            createVariableAndWrite(wds, 'theta1_tau', 'degrees',
                                'direction of first eigenvector of the symmetric part of stress gradient tensor',
                                ('time', 'Y', 'X'),
                                float,
                                np.array([theta1_tau],dtype=float))

            createVariableAndWrite(wds, 'theta2_tau', 'degrees',
                                'direction of second eigenvector of the symmetric part of stress gradient tensor',
                                ('time', 'Y', 'X'),
                                float,
                                np.array([theta2_tau],dtype=float))

            createVariableAndWrite(wds, 'eig1_vel', 'sec^-1',
                                'first eigenvalue of the symmetric part of velocity gradient tensor',
                                ('time', 'Y', 'X'),
                                float,
                                np.array([eig1_vel],dtype=float))

            createVariableAndWrite(wds, 'eig2_vel', 'sec^-1',
                                'second eigenvalue of the symmetric part of velocity gradient tensor',
                                ('time', 'Y', 'X'),
                                float,
                                np.array([eig2_vel],dtype=float))

            createVariableAndWrite(wds, 'theta1_vel', 'degrees',
                                'direction of first eigenvector of the symmetric part of velocity gradient tensor',
                                ('time', 'Y', 'X'),
                                float,
                                np.array([theta1_vel],dtype=float))

            createVariableAndWrite(wds, 'theta2_vel', 'degrees',
                                'direction of second eigenvector of the symmetric part of velocity gradient tensor',
                                ('time', 'Y', 'X'),
                                float,
                                np.array([theta2_vel],dtype=float))

            createVariableAndWrite(wds, 'curl_stress', 'dynes/cm^3',
                                'curl of stress',
                                ('time', 'Y', 'X'),
                                float,
                                np.array([omega_tau],dtype=float))

            createVariableAndWrite(wds, 'vorticity', 'sec^-1',
                                'vorticity of ocean current',
                                ('time', 'Y', 'X'),
                                float,
                                np.array([omega_vel],dtype=float))

            wds.close()

