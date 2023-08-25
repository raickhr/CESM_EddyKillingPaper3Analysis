import numpy as np
from netCDF4 import Dataset
from mpi4py import MPI
import sys


# def GetGradient(Array2D, DX, DY):
#     n = 4
#     coeffs = np.array([1./280, -4./105, 1./5, -4./5, 0, 4. /
#                       5, -1./5, 4./105, -1./280], dtype=float)

#     dx_phi = np.zeros(Array2D.shape, dtype=float)
#     dy_phi = np.zeros(Array2D.shape, dtype=float)

#     for i in range(-n, n+1):
#         dx_phi += np.roll(Array2D, i, axis=1) * -coeffs[n+i]
#         dy_phi += np.roll(Array2D, i, axis=0) * -coeffs[n+i]
#         # print(i, coeffs[n+i])

#     dx_phi = dx_phi/DX
#     dy_phi = dy_phi/DY

#     return dx_phi, dy_phi

C2dict = {'10': 1.99159,
          '20': 0.61004,
          '40': 0.25030,
          '50': 0.20579,
          '60': 0.18137,
          '80': 0.15688,
          '100': 0.14546,
          '120': 0.13923,
          '140': 0.13546,
          '150': 0.13412,
          '160': 0.13302,
          '180': 0.13134,
          '200': 0.13013,
          '220': 0.12924,
          '240': 0.12857,
          '260': 0.12804,
          '280': 0.12762,
          '300': 0.12728,
          '320': 0.12701,
          '340': 0.12678,
          '360': 0.12659,
          '380': 0.12642,
          '400': 0.12628,
          '420': 0.12617,
          '440': 0.12606,
          '500': 0.12582,
          '600': 0.12557,
          '700': 0.12542,
          '800': 0.12532,
          '900': 0.12525,
          '1000': 0.12521}


C4dict = {'10': 10.88805,
          '20': 0.90840,
          '40': 0.12277,
          '50': 0.07654,
          '150': 0.02532,
          '60': 0.05589,
          '80': 0.03850,
          '100': 0.03153,
          '120': 0.02803,
          '140': 0.02602,
          '160': 0.02475,
          '180': 0.02390,
          '200': 0.02330,
          '220': 0.02286,
          '240': 0.02253,
          '260': 0.02228,
          '280': 0.02208,
          '300': 0.02191,
          '320': 0.02178,
          '340': 0.02167,
          '360': 0.02158,
          '380': 0.02150,
          '400': 0.02144,
          '420': 0.02138,
          '440': 0.02133,
          '500': 0.02122,
          '600': 0.02110,
          '700': 0.02103,
          '800': 0.02098,
          '900': 0.02095,
          '1000': 0.02093}


def GetGradient(Array2D, DXU, DYU):
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


def getFirstAndSecondGrads(Array2D, DX, DY):
    dx, dy = GetGradient(Array2D, DX, DY)
    d2x, dydx = GetGradient(dx, DX, DY)
    dxdy, d2y = GetGradient(dy, DX, DY)
    returnArr1 = np.stack((dx, dy), axis=0)
    returnArr2 = np.stack((d2x, d2y, dxdy, dydx), axis=0)
    return returnArr1, returnArr2


def getGridData(GridFile, rank):
    ds = Dataset(GridFile)
    DXU = np.array(ds.variables['DXU'])
    DYU = np.array(ds.variables['DYU'])
    KMT = np.array(ds.variables['KMT'])

    print('Grid Read Complete at rank {0:d}'.format(rank))
    sys.stdout.flush()

    return DXU, DYU, KMT


def readFile(fileName,
             u_barVarName,
             v_barVarName,
             taux_barVarName,
             tauy_barVarName):

    ds = Dataset(fileName)
    u_bar = np.array(ds.variables[u_barVarName])
    v_bar = np.array(ds.variables[v_barVarName])
    taux_bar = np.array(ds.variables[taux_barVarName])
    tauy_bar = np.array(ds.variables[tauy_barVarName])
    ds.close()

    return u_bar, v_bar, taux_bar, tauy_bar


def readConCatFile(ds,
                   u_barVarName,
                   v_barVarName,
                   taux_barVarName,
                   tauy_barVarName,
                   day):

    index = day - 1
    u_bar = np.array(ds.variables[u_barVarName][index, :, :])
    v_bar = np.array(ds.variables[v_barVarName][index, :, :])
    taux_bar = np.array(ds.variables[taux_barVarName][index, :, :])
    tauy_bar = np.array(ds.variables[tauy_barVarName][index, :, :])

    return u_bar, v_bar, taux_bar, tauy_bar


def getModelVal(taux_bar, tauy_bar, u_bar, v_bar, DX, DY, ellinkm):
    # value for 2d plane surface , maybe calculate for spherical surface, corrected
    global C2dict, C4dict
    C_2 = C2dict[str(ellinkm)]

    dx_u, dy_u = GetGradient(u_bar, DX, DY)
    dx_v, dy_v = GetGradient(v_bar, DX, DY)

    dx_taux, dy_taux = GetGradient(taux_bar, DX, DY)
    dx_tauy, dy_tauy = GetGradient(tauy_bar, DX, DY)

    # model  =  d_j tau_i d_j u_i
    model = dx_taux * dx_u + dx_tauy * dx_v + dy_taux * dy_u + dy_tauy * dy_v

    strainVel_11 = dx_u
    strainVel_12 = 1/2 * (dx_v + dy_u)
    strainVel_21 = strainVel_12
    strainVel_22 = dy_v

    rotateVel_12 = 1/2 * (dx_v - dy_u)
    rotateVel_21 = -1*rotateVel_12

    strainTau_11 = dx_taux
    strainTau_12 = 1/2 * (dx_tauy + dy_taux)
    strainTau_21 = strainTau_12
    strainTau_22 = dy_tauy

    rotateTau_12 = 1/2 * (dx_tauy - dy_taux)
    rotateTau_21 = -1*rotateTau_12

    strainPart = strainVel_11 * strainTau_11 + strainVel_12 * strainTau_12 + \
        strainVel_21 * strainTau_21 + strainVel_22 * strainTau_22

    rotatePart = rotateVel_12 * rotateTau_12 + rotateVel_21 * rotateTau_21

    model *= 1/2 * (ellinkm*1000)**2 * C_2
    strainPart *= 1/2 * (ellinkm*1000)**2 * C_2
    rotatePart *= 1/2 * (ellinkm*1000)**2 * C_2

    return model, rotatePart, strainPart


def getModelValwith2(taux_bar, tauy_bar, u_bar, v_bar, DX, DY, ellinkm):
    global C2dict, C4dict

    # 0.125 ### value for 2d plane surface , maybe calculate for spherical surface
    M2 = C2dict[str(ellinkm)]
    M4 = C4dict[str(ellinkm)]

    d = 2

    firstGrads_u, secondGrads_u = getFirstAndSecondGrads(u_bar, DX, DY)
    firstGrads_v, secondGrads_v = getFirstAndSecondGrads(v_bar, DX, DY)

    firstGrads_taux, secondGrads_taux = getFirstAndSecondGrads(
        taux_bar, DX, DY)
    firstGrads_tauy, secondGrads_tauy = getFirstAndSecondGrads(
        tauy_bar, DX, DY)

    dx_u, dy_u = firstGrads_u[0, :, :], firstGrads_u[1, :, :]
    dx_v, dy_v = firstGrads_v[0, :, :], firstGrads_v[1, :, :]

    dx_taux, dy_taux = firstGrads_taux[0, :, :], firstGrads_taux[1, :, :]
    dx_tauy, dy_tauy = firstGrads_tauy[0, :, :], firstGrads_tauy[1, :, :]

    # model  =  d_j tau_i d_j u_i
    model = dx_taux * dx_u +\
        dy_taux * dy_u +\
        dx_tauy * dx_v +\
        dy_tauy * dy_v

    strainVel_11 = 1/2 * (dx_u + dx_u)
    strainVel_12 = 1/2 * (dx_v + dy_u)
    strainVel_21 = strainVel_12
    strainVel_22 = 1/2 * (dy_v + dy_v)

    rotateVel_12 = 1/2 * (dx_v - dy_u)
    rotateVel_21 = -1*rotateVel_12

    strainTau_11 = 1/2 * (dx_taux + dx_taux)
    strainTau_12 = 1/2 * (dx_tauy + dy_taux)
    strainTau_21 = strainTau_12
    strainTau_22 = 1/2 * (dy_tauy + dy_tauy)

    rotateTau_12 = 1/2 * (dx_tauy - dy_taux)
    rotateTau_21 = -1*rotateTau_12

    strainPart = strainVel_11 * strainTau_11 + strainVel_12 * strainTau_12 + \
        strainVel_21 * strainTau_21 + strainVel_22 * strainTau_22

    rotatePart = rotateVel_12 * rotateTau_12 + rotateVel_21 * rotateTau_21

    model *= 1/2 * (ellinkm*1e5)**2 * M2
    strainPart *= 1/2 * (ellinkm*1e5)**2 * M2
    rotatePart *= 1/2 * (ellinkm*1e5)**2 * M2

    ######################## second order ###############################
    model2 = np.zeros(np.shape(u_bar), dtype=float)
    for i in range(4):
        model2 += secondGrads_u[i, :, :]*secondGrads_taux[i, :, :] +\
            secondGrads_v[i, :, :]*secondGrads_tauy[i, :, :]

    model2 *= M4/(2*d*(d+2)) * (ellinkm*1e5)**4

    ######################## second order second part ###############################
    Del_u = secondGrads_u[0, :, :] + secondGrads_u[3, :, :]
    Del_v = secondGrads_v[0, :, :] + secondGrads_v[3, :, :]
    Del_tx = secondGrads_taux[0, :, :] + secondGrads_taux[3, :, :]
    Del_ty = secondGrads_tauy[0, :, :] + secondGrads_tauy[3, :, :]

    model3 = Del_u*Del_tx + Del_v*Del_ty

    coeff = (d*M4 - (d+2)*M2**2)/(4*d**2*(d+2))

    model3 *= coeff*(ellinkm*1e5)**4

    model2ndOrder = model2+model3
    ######################## second order complete ###############################

    return model, model2ndOrder, rotatePart, strainPart


def writeNetCDFModelAndCompsFile(writeFileName, Ylen, Xlen, timeVal, modelVal, rotPartVal, strainPartVal, ellVal, timeUnits):
    ds = Dataset(writeFileName, 'w', format='NETCDF4_CLASSIC')

    ds.createDimension('Y', Ylen)
    ds.createDimension('X', Xlen)
    ds.createDimension('time', None)
    ds.createDimension('ell', 1)

    time = ds.createVariable('time', float, ('time'))
    time.units = timeUnits
    time[0:1] = timeVal

    model = ds.createVariable('NLmodel_EPCg', float,
                              ('time', 'Y', 'X'))
    model.units = 'Watts/m^2'

    rotPart = ds.createVariable(
        'NLmodel_EPCg_rot', float, ('time', 'Y', 'X'))
    rotPart.units = 'Watts/m^2'

    strainPart = ds.createVariable(
        'NLmodel_EPCg_strain', float, ('time', 'Y', 'X'))
    strainPart.units = 'Watts/m^2'

    ell = ds.createVariable('ell', float, ('ell'))
    ell.units = 'km'

    dummy = np.zeros((1, 2400, 3600), dtype=float)
    dummy[0, :, :] = modelVal[:, :]
    model[:, :, :] = dummy[:, :, :]
    del dummy

    dummy = np.zeros((1, 2400, 3600), dtype=float)
    dummy[0, :, :] = rotPartVal[:, :]
    rotPart[:, :, :] = dummy[:, :, :]
    del dummy

    dummy = np.zeros((1, 2400, 3600), dtype=float)
    dummy[0, :, :] = strainPartVal[:, :]
    strainPart[:, :, :] = dummy[:, :, :]
    del dummy

    ell[:] = ellVal

    ds.close()


def writeNetCDFModelWith2ndOrderAndCompsFile(writeFileName, Ylen, Xlen, timeVal, modelVal, model2OrderVal, rotPartVal, strainPartVal, ellVal, timeUnits):
    ds = Dataset(writeFileName, 'w', format='NETCDF4_CLASSIC')

    ds.createDimension('Y', Ylen)
    ds.createDimension('X', Xlen)
    ds.createDimension('time', None)
    ds.createDimension('ell', 1)

    time = ds.createVariable('time', float, ('time'))
    time.units = timeUnits
    time[0:1] = timeVal

    model = ds.createVariable('NLmodel_EPCg', float,
                              ('time', 'Y', 'X'))

    model.units = 'Watts/m^2'

    model2Order = ds.createVariable('NLmodel2_EPCg', float,
                                    ('time', 'Y', 'X'))
    model2Order.units = 'Watts/m^2'

    rotPart = ds.createVariable(
        'NLmodel_EPCg_rot', float, ('time', 'Y', 'X'))
    rotPart.units = 'Watts/m^2'

    strainPart = ds.createVariable(
        'NLmodel_EPCg_strain', float, ('time', 'Y', 'X'))
    strainPart.units = 'Watts/m^2'

    ell = ds.createVariable('ell', float, ('ell'))
    ell.units = 'km'

    dummy = np.zeros((1, 2400, 3600), dtype=float)
    dummy[0, :, :] = modelVal[:, :]
    model[:, :, :] = dummy[:, :, :]
    del dummy

    dummy = np.zeros((1, 2400, 3600), dtype=float)
    dummy[0, :, :] = model2OrderVal[:, :]
    model2Order[:, :, :] = dummy[:, :, :]
    del dummy

    dummy = np.zeros((1, 2400, 3600), dtype=float)
    dummy[0, :, :] = rotPartVal[:, :]
    rotPart[:, :, :] = dummy[:, :, :]
    del dummy

    dummy = np.zeros((1, 2400, 3600), dtype=float)
    dummy[0, :, :] = strainPartVal[:, :]
    strainPart[:, :, :] = dummy[:, :, :]
    del dummy

    ell[:] = ellVal

    ds.close()


def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nprocs = comm.Get_size()

    # for rank 0 read grid file
    Xlen = 3600
    Ylen = 2400
    GridFile = '../input/tripoleGridCreated.nc'
    #newGlobalGrid_tripolePOP_0.1deg.nc'

    u_barVarName = 'UVEL'
    v_barVarName = 'VVEL'
    taux_barVarName = 'SUF'
    tauy_barVarName = 'SVF'

    # ellList = np.arange(60, 450, 20)
    # ellList = np.append(np.array([50]), ellList, axis=0)
    # ellList = np.append(ellList, np.arange(500, 1050, 100), axis=0)
    ellList = [100, 180] #, 1000]

    # DX = np.empty((Ylen, Xlen), dtype=float)
    # DY = np.empty((Ylen, Xlen), dtype=float)
    # KMT = np.empty((Ylen, Xlen), dtype=float)
    #if rank == 0:
    DX, DY, KMT = getGridData(GridFile, rank)
    
    comm.Barrier()

    # comm.Bcast(DX, root=0)
    # comm.Bcast(DY, root=0)
    # comm.Bcast(KMT, root=0)

    print('bcast grid data complete at {0:d} rank'.format(rank))
    sys.stdout.flush()

    # divide work
    nfiles = 2 #365
    nfilesForMe = nfiles//nprocs
    remainder = nfiles % nprocs

    # startIndex of filelist in each processor fileIndex start from 1
    fileListInMe = [rank + 1]  # the list goes on as [1,6,11] for 5 processors

    if rank < remainder:
        nfilesForMe += 1

    for i in range(1, nfilesForMe):
        fileListInMe.append(fileListInMe[-1]+nprocs)

    print('number of files for rank {0:d} is {1:d}'.format(rank, nfilesForMe))
    sys.stdout.flush()

    comm.Barrier()

    for ell in ellList:
        ellFldLoc = '{0:d}km_nanland/'.format(ell)
        fileConCatName = ellFldLoc + \
            'filtered_{0:04d}.nc'.format(ell)

        dsRead = Dataset(fileConCatName)
        landMask = KMT < 1

        for fileNumber in fileListInMe:
            print('reading day {0:d} in proc {1:d} for ell {2:d} km'.format(
                fileNumber, rank, ell))
            sys.stdout.flush()
            # fileName = ellFldLoc + \
            #     '{0:d}.nc_{1:d}._Filtered.nc'.format(fileNumber, ell)
            writeFileName = ellFldLoc + \
                '{0:03d}_NLmodelEP_{1:04d}km.nc'.format(fileNumber, ell)

            # u_bar, v_bar, taux_bar, tauy_bar = readFile(
            #     fileName, u_barVarName, v_barVarName, taux_barVarName, tauy_barVarName)

            u_bar, v_bar, taux_bar, tauy_bar = readConCatFile(
                dsRead, u_barVarName, v_barVarName, taux_barVarName, tauy_barVarName, fileNumber)

            #u_bar[landMask] = float('nan')
            #v_bar[landMask] = float('nan')
            #taux_bar[landMask] = float('nan')
            #tauy_bar[landMask] = float('nan')

            #u_bar[abs(u_bar) > 1e5] = float('nan')
            #v_bar[abs(v_bar) > 1e5] = float('nan')
            #taux_bar[abs(taux_bar) > 1e5] = float('nan')
            #tauy_bar[abs(tauy_bar) > 1e5] = float('nan')

            # NLmodel, rotPartModel, strainPartModel = getModelVal(taux_bar, tauy_bar,
            #                       u_bar, v_bar, DX, DY, ell)

            NLmodel, NLmodel2order, rotPartModel, strainPartModel = getModelValwith2(taux_bar, tauy_bar,
                                                                                     u_bar, v_bar, DX, DY, ell)

            timeVal = fileNumber * 7.0 + 3.0
            timeUnits = "days since 0050-01-01 00:00:00"

            # writeNetCDFModelAndCompsFile(writeFileName, lat,
            #                      lon, timeVal, NLmodel, rotPartModel, strainPartModel, ell, timeUnits)

            writeNetCDFModelWith2ndOrderAndCompsFile(writeFileName, Ylen,
                                                     Xlen, timeVal, NLmodel, NLmodel2order, rotPartModel, strainPartModel, ell, timeUnits)

        dsRead.close()

    MPI.Finalize()


if __name__ == "__main__":
    main()

