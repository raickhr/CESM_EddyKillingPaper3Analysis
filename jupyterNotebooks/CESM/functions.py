import numpy as np
from numpy import sin, cos, tan, arctan, deg2rad, rad2deg
from scipy import signal, fft, interpolate
from scipy.ndimage import gaussian_filter


def getGeodesicDistFromLonLat(center_lon, center_lat, LON, LAT, Radius):
    phi2 = (LAT)
    phi1 = (center_lat)
    dlambda = (LON - center_lon)
    num = np.sqrt((cos(phi2) * sin(dlambda))**2 +(cos(phi1)* sin(phi2) - sin(phi1)*cos(phi2)* cos(dlambda))**2)
    den = (sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2) * cos(dlambda))
    dsigma = np.arctan2(num,den)
    r = Radius*dsigma
    return r


def getXYlengthsFromLonLat(center_lon, center_lat, LON, LAT, Radius):
    dlat_fromCenter = LAT - center_lat
    dlon_fromCenter = LON - center_lon
    X = Radius * cos((LAT))*dlon_fromCenter
    Y = Radius * (dlat_fromCenter)
    return X, Y


def getKernel(rinKm, ellinKm, UAREA):
    G = (0.5 - 0.5*np.tanh((rinKm - ellinKm/2)/10))
    #G = np.ones(rinKm.shape)
    #G[rinKm > ellinKm/2] = 0.0
    normalization = np.nansum(G*UAREA)
    kernel = G/normalization
    return kernel 

def getConst(rInKm, X, Y, UAREA, ellInKm):
    kernel = getKernel(rInKm, ellInKm, UAREA)
    intgX = np.nansum(X**2 * UAREA * kernel)/(ellInKm*1e3)**2 
    intgY = np.nansum(Y**2 * UAREA * kernel)/(ellInKm*1e3)**2 
    intgXY = np.nansum(X*Y * UAREA * kernel)/(ellInKm*1e3)**2 
    intgYX = np.nansum(Y*X * UAREA * kernel)/(ellInKm*1e3)**2 
    return (intgX+intgY+intgXY+intgYX)


def getKernelAtLat(ell, lat_index, LON, LAT, UAREA, Radius):
    tolerance = 100e3
    ylen, xlen = np.shape(UAREA)
    
    x_center = int(xlen//2)
    
    center_lat = LAT[lat_index, x_center]
    center_lon = LON[lat_index, x_center]
    
    dlon = LON[0,1] - LON[0,0]
    dlat = LAT[1,0] - LAT[0,0]
    
    half_lonRange = (ell+tolerance)/(2*Radius*cos(center_lat))
    half_latRange = (ell+tolerance)/(2*Radius) 
    
    half_nlon = int(half_lonRange//dlon)
    half_nlat = int(half_latRange//dlat)
    
    x_s = x_center-half_nlon
    x_e = x_center+half_nlon
    
    if x_s < 0 or x_e > xlen:
        print('domain not wide enough !!! increase in x-direction' )
        return False
        
    y_s = lat_index - half_nlat
    y_e = lat_index + half_nlat
    
    k_center_x = half_nlon
    k_center_y = half_nlat
    
    if y_s < 0:
        y_s = 0                  ## this clips the kernel at northern border
        k_center_y = lat_index   ## y kernel center according to clipping
    if y_e > ylen:
        y_e = ylen               ## this clips the kernel at southern border
         
        
    k_UAREA = UAREA[y_s:y_e, x_s:x_e].copy()
    k_LAT = LAT[y_s:y_e, x_s:x_e].copy()
    k_LON = LON[y_s:y_e, x_s:x_e].copy()
    
    r = getGeodesicDistFromLonLat(center_lon, center_lat, k_LON, k_LAT, Radius)
    kernel = getKernel(r/1e3, ell/1e3, k_UAREA)
    
    return kernel, k_center_x, k_center_y, k_UAREA

def getFilteredAtij(field, k_UAREA, yindx, xindx, kernel, k_center_x, k_center_y):
    k_ylen, k_xlen = np.shape(kernel)
    shifted_field = np.roll(field, (k_center_y-yindx, k_center_x-xindx), axis=(0,1))
    
    k_field = shifted_field[0:k_ylen, 0:k_xlen].copy()
    field_bar_ji = np.sum(k_field * kernel * k_UAREA)
    
    return field_bar_ji


def getFilteredField(ell, field_list, UAREA, LON, LAT, Radius):
    nfields, ylen, xlen = np.shape(field_list)
    fieldbar_list = np.zeros((nfields, ylen, xlen), dtype=float)
    for j in range(ylen):
        kernel, k_center_x, k_center_y, k_UAREA = getKernelAtLat(ell, j, LON, LAT, UAREA, Radius)
        for i in range(xlen):
            for n in range(nfields):
                fieldbar_list[n,j,i] = getFilteredAtij(field_list[n,:,:], 
                                                       k_UAREA, 
                                                       j, i,
                                                       kernel, 
                                                       k_center_x, k_center_y)                  
    return fieldbar_list

    
    
    
    
def getFourierCoeffs(powerlaw, nx, ny, dx, dy):
    kx = fft.fftfreq(nx, dx)
    ky = fft.fftfreq(ny, dy)

    kX, kY = np.meshgrid(kx,ky)
    k_rad = np.sqrt(kX**2 + kY**2)
    k_max = np.max(k_rad)

    shifted_krad = fft.fftshift(k_rad)

    amp = np.empty_like(k_rad)
    amp[k_rad!=0] = k_rad[k_rad!=0]**(-powerlaw)

    amp[k_rad==0] = np.max(amp)

    coeffs = np.sqrt(amp)
    return coeffs

def getNoiseFromFourierCoeffs(coeffs, nx, ny):
    rand = np.random.uniform(low=-1,0, high=1.0, size=(ny, nx))
    #np.random.rand(ny,nx)

    randCoeffs = fft.ifft2(rand, norm='forward') 
    coeffs = coeffs*randCoeffs

    noise = fft.ifft2(coeffs, norm='forward')
    noise = noise.real
    #noise = noise/np.mean(noise)

    noise = gaussian_filter(noise, sigma=30, mode='wrap')
    noise = gaussian_filter(noise, sigma=10, mode='wrap')
    noise = gaussian_filter(noise, sigma=5, mode='wrap')
    return noise

def getNoise(powerlaw, nx, ny, dx, dy):
    coeffs = getFourierCoeffs(powerlaw, nx, ny, dx, dy)
    noise = getNoiseFromFourierCoeffs(coeffs, nx, ny)
    return noise
    
    
    

def getDXUDYU(LON, LAT, dlon_inDeg, dlat_inDeg, Radius):
    DYU = np.full(LAT.shape, Radius*deg2rad(dlat_inDeg))
    phi2 = LAT
    phi1 = LAT
    dlambda = deg2rad(dlon_inDeg)
    num = np.sqrt((cos(phi2) * sin(dlambda))**2 +(cos(phi1)* sin(phi2) - sin(phi1)*cos(phi2)* cos(dlambda))**2)
    den = (sin(phi1)*sin(phi2) + cos(phi1)*cos(phi2) * cos(dlambda))
    dsigma = np.arctan2(num,den)
    DXU = Radius*dsigma
    return DXU, DYU
    
    
    
    
    
    
    
    
    