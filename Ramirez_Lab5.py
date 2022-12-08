import glob
import rasterio
import numpy as np
from math import pi
from scipy import ndimage
import os
import scipy


# creating a file path for raster files
file_path = "C:\\Users\lreyr\Downloads\data5\data"
raster_tif = glob.glob(file_path + '/*tif')     


# opening raster as an array
elk_dem = rasterio.open(os.path.join(file_path, 'bigElk_dem.tif'))
elk_array = elk_dem.read(1)
cell_size = elk_dem.transform[0]



def slopeAspect(elk_array, cell_size):
    """Calculates slope and aspect using the 3rd-order finite difference method
    Parameters
    ----------
    dem : numpy array
        A numpy array of a DEM
    cs : float
        The cell size of the original DEM
    Returns
    -------
    numpy arrays
        Slope and Aspect arrays
    """
    from math import pi
    from scipy import ndimage
    kernel = np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]])
    dzdx = ndimage.convolve(elk_array, kernel, mode='mirror') / (8 * cell_size)
    dzdy = ndimage.convolve(elk_array, kernel.T, mode='mirror') / (8 * cell_size)
    slp = np.arctan((dzdx ** 2 + dzdy ** 2) ** 0.5) * 180 / pi
    ang = np.arctan2(-dzdy, dzdx) * 180 / pi
    aspect = np.where(ang > 90, 450 - ang, 90 - ang)
    return slp, aspect


elk_slope, elk_aspect = slopeAspect(elk_array, cell_size)

def reclassAspect(npArray):
    """Reclassify aspect array to 8 cardinal directions (N,NE,E,SE,S,SW,W,NW),
    encoded 1 to 8, respectively (same as ArcGIS aspect classes).
    Parameters
    ----------
    npArray : numpy array
        numpy array with aspect values 0 to 360
    Returns
    -------
    numpy array
        numpy array with cardinal directions
    """
    return np.where((npArray > 22.5) & (npArray <= 67.5), 2,
    np.where((npArray > 67.5) & (npArray <= 112.5), 3,
    np.where((npArray > 112.5) & (npArray <= 157.5), 4,
    np.where((npArray > 157.5) & (npArray <= 202.5), 5,
    np.where((npArray > 202.5) & (npArray <= 247.5), 6,
    np.where((npArray > 247.5) & (npArray <= 292.5), 7,
    np.where((npArray > 292.5) & (npArray <= 337.5), 8, 1)))))))

# ReclassAspect using aspect function
reclassAspect(elk_slope)


def reclassByHisto(npArray, bins):
    """Reclassify np array based on a histogram approach using a specified
    number of bins. Returns the reclassified numpy array and the classes from
    the histogram.
    Parameters
    ----------
    npArray : numpy array
        Array to be reclassified
    bins : int
        Number of bins
    Returns
    -------
    numpy array
        umpy array with reclassified values
    """
    # array = np.where(np.isnan(npArray), 0, npArray)
    histo = np.histogram(~np.isnan(npArray), bins)[1]
    rClss = np.zeros_like(npArray)
    for i in range(bins):
        print(i + 1, histo[i], histo[i + 1])
        print(np.where((npArray > histo[i]) & (npArray <= histo[i + 1])))
        rClss = np.where((npArray >= histo[i]) & (npArray <= histo[i + 1]),
                         i + 1, rClss)
    return rClss


# Reclassify np array using reclassByHisto function
reclassByHisto(elk_slope, 10)

# Open fire raster to use during recovery ratio by year
fire_rast = rasterio.open(os.path.join(file_path, 'fire_perimeter.tif'))
fire_array = fire_rast.read(1)
fire_profile = fire_rast.profile



# create a list for band3 and band4 files
from pathlib import Path

path_of_the_directory = "C:\\Users\lreyr\Downloads\data5\data\L5_big_elk"
file = Path(path_of_the_directory ).glob('*B3.tif')
for B3 in file:
    print(B3)

path_of_the_directory2 = "C:\\Users\lreyr\Downloads\data5\data\L5_big_elk"
file = Path(path_of_the_directory ).glob('*B4.tif')
for B4 in file:
    print(B4)


# Open and create arrays for b3 and b4 files
band3 = rasterio.open(B3)
band4 = rasterio.open(B4)
band3_array = band3.read(1)
band4_array = band4.read(1)

# calculating NDVI for images with bands 3 and 4
NDVI = (band4_array - band3_array) / (band4_array + band3_array)


# calculating recovery ratio for each pizel for each year
meantemp = np.mean(NDVI[np.where(fire_array == 2)])
recovery_ratio = NDVI / meantemp








