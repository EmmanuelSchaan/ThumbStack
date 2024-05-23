import numpy as np
from tqdm import tqdm

def f(xR, cTr): 
    """This function appears in the 
    lensing deflection angle profile
    of an NFW halo density truncated at Rvir
    xR & cTr: both dimensionless"""
    
    mask1 = xR.real < cTr
    mask2 = xR.real > 0.01
    mask3 = xR.real < 0.01
    
    f = np.zeros(xR.shape, dtype=np.complex128)
    
    f[mask1 & mask2] = -xR[mask1 & mask2]/(1+cTr)/(cTr+np.sqrt(cTr**2-xR[mask1 & mask2]**2)) \
        +(np.log(1+cTr)+np.log(xR[mask1 & mask2])-np.log(cTr+np.sqrt(cTr**2-xR[mask1 & mask2]**2)))/xR[mask1 & mask2] \
        -np.log(-cTr-xR[mask1 & mask2]**2+np.sqrt(1-xR[mask1 & mask2]**2)*np.sqrt(cTr**2-xR[mask1 & mask2]**2))/xR[mask1 & mask2]/np.sqrt(1-xR[mask1 & mask2]**2) \
        +(np.log(1+cTr)+np.log(-xR[mask1 & mask2])-2j*(np.log(-xR[mask1 & mask2]).imag))/xR[mask1 & mask2]/np.sqrt(1-xR[mask1 & mask2]**2)
    
    f[mask1 & mask3] = 0 + 0j
    
    f[~mask1 & mask2] = (1/xR[~mask1 & mask2]) * ((-cTr/(1+cTr)) + np.log(1+cTr))
        
    return f.real


def paintML(datFram):
    """
    paints the moving lens temperature anisotropy
    on a flat patch of sky given the catalog of objects
    """
    
    "flat map"
    # number of pixels for the flat map
    nX = 3*2048 #3414 #4779 #3*2048 #4096
    # nY = 1* 2048 #4096
    # flat map and catalog patch dimensions in degrees
    sizeX = 45 *np.pi/180.
    # sizeY = 45
    # make a plain flat map
    dX = float(sizeX)/(nX-1)
    x1d = dX * np.arange(nX)   # the x value corresponds to the center of the cell
    # dY = float(sizeY)/(nY-1)
    # y = dY * np.arange(nY)   # the y value corresponds to the center of the cell
    x0, y0 = np.meshgrid(x1d, x1d, indexing='ij')
    # get the resolution of flat map in arcmin
    baseMapRes = dX *60*180/np.pi #arcmin
        
    "general constants"
    c_ms = 299792458.0   # Speed of light in m/s
    G = 0.004299421976550123 # [m^2.Mpc.c^-2.Msun^-1]
    c_kms = 299792.458  # Speed of light in km/s
    Tcmb= 2.726e6  # uK

    "halos info"
    z = datFram.Z # redshift
    a = 1/(1+z) # scale factor
    cNFW = datFram.cNFW 
    Rs = datFram.Rs # (comoving) [Mpc]
    rhoS = datFram.rhoS # (comoving) [Msun/Mpc^3]
    cTr = cNFW # truncation radius in units of Rs
    comovDist = datFram.comovDist # (comoving) Mpc
    vTh = datFram.vTh # km/s
    vPh = datFram.vPh # km/s
    RArad = datFram.RArad
    DECrad = datFram.DECrad

    "canvas"
    # make a plain canvas
    canvas = np.zeros(np.shape(x0))

    # "spline"
    # # sample
    # xsamp = np.linspace(0.1, 14670, 100)
    # xsamp = xsamp.astype(complex)
    # ysamp = f(xsamp, 4.826226)
    # # spline
    # xsamp = np.linspace(0.1, 14670, 100)
    # spline_fnf = interp1d(xsamp, ysamp, kind='cubic')
    # "adaptive sampled spline"
    # # Initial sampling
    # initial_samples = np.linspace(0.1, 60864, 1000)
    # function_values = f(initial_samples.astype(complex), 4.826226)
    # # Identify regions of interest (for example, where the function changes rapidly)
    # variation_indices = np.where(np.abs(np.diff(function_values)) > 0.1)[0]
    # # Define a function to calculate the density of additional samples
    # def sample_density(x):
    #     return 50.0 + np.exp(-x)  # Adjust the function as needed
    # # Adaptive sampling in regions of interest
    # additional_samples = []
    # for idx in variation_indices:
    #     start = initial_samples[idx]
    #     end = initial_samples[idx + 1]
    #     density = int(sample_density((end - start) / 2))
    #     new_samples = np.linspace(start, end, density + 1)
    #     additional_samples.extend(new_samples)
    # # Combine initial and additional samples
    # all_samples = np.sort(np.concatenate((initial_samples, additional_samples)))
    # # Evaluate the function and construct splines
    # interp_function = interp1d(all_samples, f(all_samples.astype(complex),4.826226))

    "draw ML map for each object and spray on canvas"
    for i in tqdm(range(len(datFram.Z))):
        
        # get x & y projections of flat map
        x = np.copy(x0) #rad
        y = np.copy(y0) #rad

            # step3) move the origin of x & y maps to the halo's position
        x -= DECrad[i] #rad
        y -= RArad[i] #rad

        # make the radius and angle maps [centered at halo's position]
        r = np.sqrt(y**2 + x**2) #rad
        angle = np.arctan2(x,y) #rad

        # convert the radius map to (comoving)Mpc
        rComov = r *comovDist[i] #(comoving)Mpc

        # make the scaled radius map [dimensionless]
        xR = rComov /Rs[i] #dimensionless
        xR = xR.astype(complex)

        # map out the f on scaled radius map
        notorF = f(xR, 1*cTr[i])

        # make the deflection angle (beta) map
        beta = ( 16. *np.pi *G *rhoS[i] *(Rs[i]**2) /(a[i]*(c_ms**2)) ) *(notorF) # dimensionless

        # find the direction of halo's transverse velocity 
        vAngle = np.arctan2(vTh[i],vPh[i])

        # make the moving lens dT map
        dT = beta * np.sqrt(vTh[i]**2 +vPh[i]**2) * np.cos(vAngle-angle) /c_kms *Tcmb # uK

        # spray dT map of each halo onto the canvas
        canvas += dT

    print("map resolution:", dX*60*180/np.pi)
    return canvas