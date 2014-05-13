from __future__ import division

import numpy as np
import scipy.ndimage as ndimage
from astropy.io import fits

from misc import *


def extract_island(coord, data, n=1):
    '''Extract an island from an array

    This function extract a small island of data from a numpy array. It
    deals with edge effects, such that resulting array always has the same
    size and shape. The result will be quadratic with dimension [2n+1, 2n+1]
    and the central pixel corresponds to the coordinates requested.
    
    Parameters
    ----------
    coord : tuple or array
        x, y coordinates (in python convention, i.e. (0,0) is the first pixel)
    data : np.array
        data array, can be masked.
    n : integer
        size of the output array is [2n+1, 2n+1]

    Returns
    -------
    island : np.ma.array
        Values that are masked in input or that go beyond the input data
        on the edges are masked.
    '''
    x, y = coord
    island = np.ma.zeros((2*n+1, 2*n+1))
    island[:] = np.ma.masked
    x0 = max(0, x-n)
    x1 = min(data.shape[0], x+n+1)
    y0 = max(0, y-n)
    y1 = min(data.shape[1], y+n+1)

    island[max(0,n-x):min(2*n+1,data.shape[0]-x+n),
           max(0,n-y):min(2*n+1,data.shape[1]-y+n)] = data[x0:x1, y0:y1]
    return island


def write_island(coord, data, island):
    '''Write to a small island in an array

    This function replaces the values in  a small island of data in a numpy 
    array. It deals with edge such that data beyond the range of the numpy
    array is ignored.
    
    Parameters
    ----------
    coord : tuple or array
        x, y coordinates (in python convention, i.e. (0,0) is the first pixel)
    data : np.array
        data array, can be masked.
    island : np.ma.array
        Values that are masked in input or that go beyond the input data
        on the edges are masked.

    Returns
    -------
    data : array
        modified from input
    '''

    if not np.all(np.array(island.shape) == island.shape[0]):
        raise ValueError('island must have quadratic shape')
    if not np.mod(island.shape[0],2) == 1:
        raise ValueError('island must have odd number elements per dimension')
    x, y = coord
    n = island.shape[0]//2-1
    x0 = max(0, x-n)
    x1 = min(data.shape[0], x+n+1)
    y0 = max(0, y-n)
    y1 = min(data.shape[1], y+n+1)

    data[x0:x1, y0:y1] = island[max(0,n-x):min(2*n+1,data.shape[0]-x+n),
           max(0,n-y):min(2*n+1,data.shape[1]-y+n)]
    return data



class BaseReplacer(object):
    '''An object that controls replacing bad pixel values in a datafile

    The BaseReplacer itself only provides core methods that apply to each
    type of replacement. It needs to be extended with a specific method
    to determine the new values that shall overwrite the bad pixels.
    '''
    replace_mask = None

    def __call__(self, data, error, replace_mask):
        '''error can be None e.g. in raw files'''
        self.data = data
        self.error = error
        self.replace_mask = replace_mask
        self.replace_all()
        return self.data, self.error, self.replace_mask


    def run_on_hdu(self, hdus, oref='.', hotthreshold=1e-4, iter=1):
        '''Run replacement on a hdulist

        Parameters
        ----------
        hdus : list of fits.hdu objects
            This is the raw file read by astropy.fits.open
        oref : string
            The path to certain calibration files (e.g. the dark) is given
            in the header of the hdus in IRAF notation. oref is the variable
            that contains the path to these files. 
        hotthreshold : float
            Pixels with a flux value above hotthreshold will be masked for
            replacement.
        iter : integer
            Number of times the replacer is run in a loop. This can be useful,
            if data in an 3*3 island is masked a invalid. In the first pass, 
            the eight outside pixels with be replaced, in the second pass, the
            central pixel is set to the mean of its eight neighbors.
        ''' 
        data = hdus[1].data
        error = hdus[2].data
        darkfile = IRAF2filepath(hdus[0].header['DARKFILE'], oref=oref)
        mask = fits.getdata(darkfile, 1) > hotthreshold
        for i in range(iter):
            data, error, mask = self(data, error, mask)
        if self.replace_mask.sum() > 0:
            print '{0} masked pixels remain uncorrected'.format(mask.sum())
        hdus[1].data = data
        if error is not None:
            hdus[2].data = error
     

    def replace_xy(self, x,y):
        raise NotImplementedError


    def replace_all(self):
        '''run a loop over all bad pixels and call replace_xy for each of them'''
        if self.replace_mask is not None:
            for x,y in zip(*self.replace_mask.nonzero()):
                replace_data, replace_err = self.replace_xy(x,y)
                if replace_data is not None:
                    self.data[x,y] = replace_data
                    if self.error is not None:
                        self.error[x,y] = replace_err
                    self.replace_mask[x,y] = False
        else:
            print "Set object.replace_mask to specify points to be replaced."


class SurroundingMeanReplacer(BaseReplacer):
    def replace_xy(self, x, y):
       '''replace bad pixels with the mean of the valid neighboring pixels'''
        dataisland = extract_island((x, y), 
                           np.ma.array(self.data, mask=self.replace_mask))
        datamean = dataisland.mean()
        if self.error is not None:
            errorisland = extract_island((x, y),
                           np.ma.array(self.error, mask=self.replace_mask))
            errormean = errorisland.mean()
        else:
            errormean = None
        return datamean, errormean


class BaseSmoother(object):
    '''Smoothes data
    '''
    def __call__(self, data, **kwargs):
        return ndimage.gaussian_filter(data, **kwargs)

    def run_on_hdu(self, hdus):
        hdus[1].data = self(hdus[1].data)
        

class InterOrderSmoother(BaseSmoother):
    '''select the regions between spectral orders and smooth only the
    data in those in a flux conserving manner. This code is in progress and not
    usable yet:
    - the make_interordermask does not select the right regions (probably that
      has something to do with the fact that the processed images where I 
      read the keywords are rebinned compared to the raw files).
    - the smoothing is FAR to slow.
    '''
    def run_on_hdu(self, hdus, hdus_x1d, n, sigma):
        self.make_interordermask(hdus, hdus_x1d)
        data = self(hdus[1].data, self.interordermask, n, sigma)
        hdus[1].data = data

    def make_interordermask(self, hdus, hdus_x1d):
        '''make a mask that covers all regions that contain the object spectrum

        Currently, this does not work

        Parameters
        ----------
        hdus : list of fits.hdu objects
            This is the raw file read by astropy.fits.open
        hdus_x1d : list of fits.hdu objects
            This is the x1d file read by astropy.fits.open
            In the x1d certain header keywords contain the position of the 
            spectral orders. Instead of rerunning the fitting tasks for that
            we just read the results from that pipeline processed file.

        '''
        data = hdus_x1d[1].data
        self.interordermask = np.zeros_like(hdus[1].data, dtype=np.bool)
        SHIFTA2 = hdus_x1d[1].header['SHIFTA2']
        CRSCROFF = hdus_x1d[1].header['CRSCROFF']
        for i in range(len(data)):
            EXTRSIZE = data[i]['EXTRSIZE']
            y0 = data[i]['A2CENTER']+SHIFTA2+CRSCROFF+EXTRSIZE/2.
            # +1 because fits convention is to start with (1,1), while numpy
            # uses (0,0).
            self.interordermask[2*np.floor(y0)+1:2*np.ceil(y0+EXTRSIZE)+1,:] = True

    def kernel(self, sigma, n):
        """ Returns a normalized 2D gauss kernel array for convolutions 

        This kernel is normalized such that the sum over all elements is 
        one, which means that an integral on the infinite plane would be > 1.
        """
        x, y = np.mgrid[-n:n+1, -n:n+1]
        g = np.exp(-(x**2+y**2)/(2.*sigma**2.))
        return g / g.sum()
        
    def __call__(self, data, mask, n, sigma):
        mdata = np.ma.array(data, mask=mask, copy=True) 
        # set all bg points to zero and build up again
        data[mask] = 0.
        kernel = self.kernel(sigma, n)
        for x,y in zip(*mask.nonzero()):
            island = extract_island((x,y), mdata, n=n)
            kmask = np.ma.array(kernel, mask=island.mask)
            kmask = kmask / kmask.sum()
            data += write_island((x,y), np.zeros_like(data), kmask*data[x,y])
        return data

