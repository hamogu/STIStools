
cp *raw.fits no_bkg/
cp *wav.fits no_bkg/

python
import glob

from astropy.io import fits

flist = glob.glob('*raw.fits')
for f in flist:
    hdulist = fits.open(f, mode='update')
    hdulist[0].header['BACKCORR'] = 'OMIT    '
    hdulist[0].header['SC2DCORR'] = 'OMIT    '
    hdulist.close()
import os

class BaseReplacer(object):
    replace_mask = None
    def replaceval_xy(self, x,y):
        raise NotImplementedError
    def set_replace_mask(self, mask):
        self.replace_mask = mask
    def replace_all(self):
        if self.replace_mask:
            for x,y in self.replace_mask.nonzero():
                replace_data, replace_err = self.replaceval_xy(x,y)
                if replace_val:
                    self.data[x,y] = replace_data
                    self.error[x,y] = replace_err
                    self.replace_mask[x,y] = False
        else:
            print "Use set_replace_mask to specify points to be replaced."

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
    if not data.shape==mask.shape:
        raise ValueError('Data and mask must have same shape.')
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




class SurroundgingMeanReplacer(BaseReplacer):
    def replace_xy(self, x, y):
        dataisland = extract_island(x, y, self.data, self.replace_mask)
        errorisland = extract_island(x, y, self.error, self.replace_mask) 
        return dataisland.mean(), errorisland.mean()



def replace_bad_pixels(self, hdus, oref='.', hotthreshold=1e-4):
    data = hdus[1].data
    darkfile = IRAF2filepath(hdus[0].header['DARKFILE'])
    mask = darkfile[1].data > hotthreshold
    for x, y in mask.nonzero():
        data[x,y] = self.Replacer(data, x,y)
    
